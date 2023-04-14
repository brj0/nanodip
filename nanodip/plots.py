"""
## Plots

Functions for creating Copy Number Variation plot.
Functions for creating methylation UMAP plot.
"""

# start_external_modules
import bisect
import csv
import logging
import os

from plotly.io import write_json, from_json
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
# end_external_modules

# start_internal_modules
from nanodip.config import (
    CNV_GRID,
    CNV_LINK,
    ENDING,
    NANODIP_REPORTS,
    PLOTLY_RENDER_MODE,
    UMAP_PLOT_TOP_MATCHES,
)
from nanodip.data import (
    Reference,
    Genome,
    Sample,
    get_reference_methylation,
    get_sample_methylation,
)
from nanodip.utils import (
    bonferroni_corrected_ci,
    convert_html_to_pdf,
    discrete_colors,
    render_template,
    composite_path,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)

def get_bin_edges(n_bins, genome):
    """Returns sequence of {n_bins} equal sized bins on chromosomes used
    by numpy histogram. Every bin is limited to one chromosome.
    """
    edges = np.linspace(0, len(genome), num=n_bins + 1).astype(int)
    # limit bins to only one chromosome
    for chrom_edge in genome.chrom.offset:
        i_nearest = np.abs(edges - chrom_edge).argmin()
        edges[i_nearest] = chrom_edge
    return edges

def get_cnv(read_positions, genome):
    """Returns CNV.

    Args:
        read_positions: List of reads in the form [start, end].
        genome: Reference genome.

    Returns:
        bin_midpoints: numpy array of x-values
        copy_numbers: numpy array of y-values
    """
    expected_reads_per_bin = 30
    n_bins = len(read_positions)//expected_reads_per_bin
    read_start_positions = [i[0] for i in read_positions]
    copy_numbers, bin_edges = np.histogram(
        read_start_positions,
        bins=get_bin_edges(n_bins, genome),
        range=[0, len(genome)],
    )
    bin_midpoints = (bin_edges[1:] + bin_edges[:-1])/2
    return bin_midpoints, copy_numbers

def cnv_grid(genome):
    """Returns chromosome grid layout for CNV Plot as plotly object and
    saves it on disk. If available grid is directly read from disk.
    """
    # Check if grid exists and return if available.
    if os.path.exists(CNV_GRID):
        with open(CNV_GRID, "r") as f:
            grid = from_json(f.read())
        return grid

    grid = go.Figure()
    grid.update_layout(
        coloraxis_showscale=False,
        xaxis = dict(
            linecolor="black",
            linewidth=1,
            mirror=True,
            range=[0, len(genome)],
            showgrid=False,
            ticklen=10,
            tickmode="array",
            ticks="outside",
            tickson="boundaries",
            ticktext=genome.chrom.name,
            tickvals=genome.chrom.center,
            zeroline=False,
        ),
        yaxis = dict(
            linecolor="black",
            linewidth=1,
            mirror=True,
            showline=True,
        ),
        template="simple_white",
    )
    # Vertical line: centromere.
    for i in genome.chrom.centromere_offset:
        grid.add_vline(x=i, line_color="black", line_dash="dot", line_width=1)
    # Vertical line: chromosomes.
    for i in genome.chrom.offset.tolist() + [len(genome)]:
        grid.add_vline(x=i, line_color="black", line_width=1)
    # Save to disk
    grid.write_json(CNV_GRID)
    return grid

def cnv_plot_from_data(data_x, data_y, expt_y, sample_name, read_num, genome):
    """Create CNV plot from CNV data.

    Args:
        data_x: x-Values to plot.
        data_y: y-Values to plot.
        expt_y: expected y-Value.
        sample_name: Name of sample.
        read_num: Number of read reads.
        genome: Reference Genome.
    """
    grid = cnv_grid(genome)
    # Expected value: draw horizontal line.
    grid.add_hline(y=expt_y, line_color="black", line_width=1)
    plot = px.scatter(
        x=data_x,
        y=data_y,
        labels={
            "x":f"Number of mapped reads: {read_num}",
            "y":f"Copy numbers per {round(len(genome)/(len(data_x)*1e6), 2)} MB"
        },
        title=f"Sample ID: {sample_name}",
        color=data_y,
        range_color=[expt_y*0, expt_y*2],
        color_continuous_scale="Portland",
        render_mode=PLOTLY_RENDER_MODE,
    )
    plot.update_traces(hovertemplate="Copy Numbers = %{y} <br>")
    plot.update_layout(grid.layout, yaxis_range = [-0.5, 2*expt_y])
    return plot

def number_of_reads(sorted_read_start_pos, interval):
    """Return the number of starting sequences within interval. Reads must
    be sorted in ascending order.
    """
    left, right = interval
    i_left = bisect.bisect_left(sorted_read_start_pos, left)
    i_right = bisect.bisect_left(sorted_read_start_pos, right)
    return len(sorted_read_start_pos[i_left:i_right])

def cnv_plot(sample, bin_midpoints, cnv, genome):
    """Create a genome-wide copy number plot and save data on dist."""
    logger.info("CNVP start")
    logger.info(sample)
    logger.info("Bin midpoints:\n%s", bin_midpoints)
    logger.info("CNV:\n%s", cnv)

    avg_read_per_bin = len(sample.reads) // len(bin_midpoints)

    plot = cnv_plot_from_data(
        data_x=bin_midpoints,
        data_y=cnv,
        expt_y = avg_read_per_bin,
        sample_name=sample.name,
        read_num=len(sample.reads),
        genome=genome,
    )
    logger.info("CNVP done")
    return plot

class CNVData:
    """CNV data container and methods for invoking CNV plot algorithm."""
    genome = Genome()
    def __init__(self, sample_name):
        self.sample = Sample(sample_name)
        self.plot = None
        self.plot_json = None
        self.genes = None
        self.relevant_genes = None
        self.bin_midpoints = None
        self.cnv = None

    def path(self, ending):
        """Returns generic path with corresponding ending."""
        return composite_path(
            NANODIP_REPORTS, self.sample.name, ENDING[ending]
        )

    def files_on_disk(self):
        """Checks if files are on disk."""
        return (
            os.path.exists(self.path("binmdpnts_npy")) and
            os.path.exists(self.path("cnv_npy")) and
            os.path.exists(self.path("cnv_json")) and
            os.path.exists(self.path("genes_csv"))
        )

    def read_from_disk(self):
        """Reads files from disk."""
        with open(self.path("cnv_json"), "r") as f:
            self.plot_json = f.read()
        self.plot = from_json(self.plot_json)
        self.genes = pd.read_csv(self.path("genes_csv"))
        self.bin_midpoints = np.load(
            self.path("binmdpnts_npy"), allow_pickle=True,
        )
        self.cnv = np.load(
            self.path("cnv_npy"), allow_pickle=True,
        )

    def make_cnv_plot(self):
        """Generates CNV plot and saves to disk."""
        self.sample.set_reads() # time consumption 9s
        self.bin_midpoints, self.cnv = get_cnv(
            self.sample.reads,
            CNVData.genome,
        )
        if len(self.bin_midpoints) == 0:
            raise ValueError("no points to plot")
        self.plot = cnv_plot(
            sample=self.sample,
            bin_midpoints=self.bin_midpoints,
            cnv=self.cnv,
            genome=CNVData.genome,
        )
        self.plot_json = self.plot.to_json()
        self.genes = self.gene_cnv()
        self.relevant_genes = self.genes.loc[self.genes.relevant]
        self.save_to_disk()

    def save_to_disk(self):
        """Saves relevant data to disk."""
        np.save(self.path("binmdpnts_npy"), self.bin_midpoints)
        np.save(self.path("cnv_npy"), self.cnv)
        self.plot.write_html(
            self.path("cnv_html"),
            config=dict({"scrollZoom": True}),
        )
        write_json(self.plot, self.path("cnv_json"))
        # time consuming operation (11s)
        self.plot.write_image(
            self.path("cnv_png"), width=1280, height=720, scale=3,
        )
        with open(self.path("alignedreads_txt"), "w") as f:
            f.write(f"{len(self.sample.reads)}")
        with open(self.path("reads_csv"), "w") as f:
            write = csv.writer(f)
            write.writerows(self.sample.reads)
        self.genes.to_csv(self.path("genes_csv"), index=False)
        self.relevant_genes.to_csv(self.path("relgenes_csv"), index=False)

    def gene_cnv(self):
        """Returns pandas DataFrame containing copy number variation
        for all genes in reference genome.
        """
        genes = CNVData.genome.genes
        genes["interval"] = list(zip(genes.start, genes.end))
        read_start_pos = [i[0] for i in self.sample.reads]
        read_start_pos.sort()
        genes["cn_obs"] = genes.interval.apply(
            lambda z: number_of_reads(read_start_pos, z)
        )
        genes.drop("interval", axis=1, inplace=True)
        bin_size = len(CNVData.genome)/(len(self.bin_midpoints))
        genes["cn_per_bin"] = genes.apply(
            lambda z: z["cn_obs"]/z["len"] * bin_size, # TODO auto draw extreme values
            axis=1,
        )
        genes["ci_left"], genes["ci_right"] = bonferroni_corrected_ci(
            hits=genes["cn_obs"],
            lengths=genes["len"],
            trials=len(read_start_pos),
            target_length=bin_size,
        )
        genes["cn_exp"] = genes.apply(
            lambda z: len(self.sample.reads)*z["len"]/len(CNVData.genome),
            axis=1,
        )
        genes["cn_obs_exp_ratio"] = genes.apply(
            lambda z: z["cn_obs"]/z["cn_exp"],
            axis=1,
        )
        genes = genes.sort_values(by="cn_obs", ascending=False)
        return genes

    def get_gene_positions(self, genes):
        """Returns sub-DataFrame for the genes of the list {genes}.
        """
        gene_pos = self.genes.loc[self.genes.name.isin(genes)]
        return gene_pos

    def plot_cnv_and_genes(self, gene_names=[]):
        """Returns json plot of the CNV plot including the CN of all
        genes in the list {gene_names}.
        """
        genes = self.get_gene_positions(gene_names)
        bin_size = len(CNVData.genome)/(len(self.bin_midpoints))
        plot = go.Figure(self.plot)
        plot.add_trace(
            go.Scatter(
                customdata=genes[[
                    "name",          # 0
                    "loc",           # 1
                    "transcript",    # 2
                    "len",           # 3
                    "cn_obs",        # 4
                    "cn_exp",        # 5
                ]],
                hovertemplate=(
                    "<b> %{customdata[0]} </b> "
                    "%{customdata[1]} "
                    "(hg19 %{customdata[2]}) <br>"
                    "Copy numbers per "
                    f"{round(bin_size/1e6, 2)}"
                    " MB: %{y} <br>"
                    "Gene length: %{customdata[3]} base pairs <br>"
                    "Observed hits: %{customdata[4]}<br>"
                    "Expected hits: %{customdata[5]:.3f}<br>"
                ),
                name="",
                marker_color="rgba(0,0,0,1)",
                mode="markers+text",
                marker_symbol="diamond",
                textfont_color="rgba(0,0,0,1)",
                showlegend=False,
                text=genes.name,
                textposition="top center",
                x=genes.midpoint,
                y=genes.cn_per_bin,
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=genes.ci_right - genes.cn_per_bin,
                    arrayminus=genes.cn_per_bin - genes.ci_left,
                ),
            ))
        return plot.to_json()

def umap_plot_from_data(sample, reference, umap_df, close_up):
    """Create and return umap plot from UMAP data.

    Args:
        sample: Sample data.
        reference: Reference data.
        umap_df: pandas data frame containing UMAP matrix and
            attributes. First row,w corresponds to sample.
        close_up: Bool to indicate if only top matches should be plotted.
    Returns:
        UMAP plot as plotly object.
    """
    # If true, sample methylation is part of analysis.
    add_sample = not sample.cpgs_only()

    umap_sample = umap_df.iloc[0]
    title0 = f"for {sample.name}" if add_sample else ""
    umap_title = (
        f"UMAP {title0} <br><sup>Reference: {reference.name} "
        f"({len(reference.specimens)} cases), "
        f"{len(sample.cpg_overlap)} CpGs </sup>"
    )
    if close_up:
        umap_title = "Close-up " + umap_title
    methyl_classes = umap_df.methylation_class[1:].to_list()
    methyl_classes.sort()
    umap_plot = px.scatter(
        umap_df,
        x="x",
        y="y",
        labels={"x":"UMAP 0", "y":"UMAP 1", "methylation_class":"WHO class"},
        title=umap_title,
        color="methylation_class",
        color_discrete_map={
            sample.name: "#ff0000",
            **discrete_colors(methyl_classes),
        },
        hover_name="id",
        category_orders={"methylation_class": [sample.name] + methyl_classes},
        hover_data=["description"],
        render_mode=PLOTLY_RENDER_MODE,
        template="simple_white",
    )
    if add_sample:
        umap_plot.add_annotation(
            x=umap_sample["x"],
            y=umap_sample["y"],
            text=sample.name,
            showarrow=True,
            arrowhead=1,
        )
    umap_plot.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
        mirror=True,
    )
    umap_plot.update_xaxes(
        mirror=True,
    )
    umap_df_ref = umap_df[1:] if add_sample else umap_df
    links = [
        f"<a href='{url}' target='_blank'>&nbsp;</a>"
        for url in [CNV_LINK % id_ for id_ in umap_df_ref.id]
    ]
    # Add hyperlinks
    umap_plot.add_trace(go.Scatter(
        x=umap_df_ref.x,
        y=umap_df_ref.y,
        mode="text",
        name="CNV links",
        text=links,
        hoverinfo="skip",
        visible="legendonly",
    ))
    # If close-up add hyperlinks for all references and draw circle
    if close_up:
        umap_plot.update_traces(marker=dict(size=5))
        # Draw circle
        radius = umap_df["distance"].iloc[-1]
        umap_plot.add_shape(
            type="circle",
            x0=umap_sample["x"] - radius,
            y0=umap_sample["y"] - radius,
            x1=umap_sample["x"] + radius,
            y1=umap_sample["y"] + radius,
            fillcolor="rgba(0,0,0,0)",
            line_color="black",
            line_width=1.0,
        )
    return umap_plot

def dimension_reduction(reference, sample):
    """Performs UMAP 2d-dimension reduction.

    Args:
        reference: Reference data which are compared with sample.
        sample: Sample to analyse containing CpG's possibly with
            methylation status.
    Returns:
        methyl_mtx: Methlation matrix with rows corresponding to
            reference specimens and columns corresponding to CpGs of
            sample overlapping reference. If sample contains methylation
            status for each CpG, then a first row corresponding to
            sample is added to matrix.
        umap_df: UMAP DataFrame corresponding to dimension reduction of
            methyl_mtx.
    """
    # Moved here due to long loading time (13.5s)
    import umap

    # If true, sample methylation is part of analysis.
    add_sample = not sample.cpgs_only()

    logger.info("Start UMAP for %s / %s.", sample.name, reference.name)
    logger.info(reference)

    # Calculate overlap of sample CpG's with reference CpG's (some probes
    # have been skipped from the reference set, e.g. sex chromosomes).
    logger.info(sample)

    if not sample.cpg_overlap:
        logger.info("UMAP done. No Matrix created, no overlapping data.")
        raise ValueError("Sample has no overlapping CpG's with reference.")

    # Extract reference and sample methylation according to CpG overlap.
    reference_methylation = get_reference_methylation(
        reference.specimens, sample.cpg_overlap,
    )
    logger.info("Reference methylation extracted:\n%s", reference_methylation)

    if add_sample:
        sample_methylation = get_sample_methylation(sample)
        logger.info("Sample methylation extracted:\n%s", sample_methylation)
        methyl_mtx = np.vstack([sample_methylation, reference_methylation])
    else:
        methyl_mtx = reference_methylation

    # Calculate UMAP Nx2 Matrix. Time intensive (~1min).
    logger.info("UMAP algorithm initiated.")
    umap_2d = umap.UMAP(verbose=True).fit_transform(methyl_mtx)
    logger.info("UMAP algorithm done.")

    # Free memory
    del reference_methylation

    if add_sample:
        umap_sample = umap_2d[0]
        umap_df = pd.DataFrame({
            "distance": [np.linalg.norm(z - umap_sample) for z in umap_2d],
            "methylation_class": [sample.name] + reference.methylation_class,
            "description":  ["Analysis sample"] + reference.description,
            "id": [sample.name] + reference.specimens,
            "x": umap_2d[:,0],
            "y": umap_2d[:,1],
        })
    else:
        umap_df = pd.DataFrame({
            "methylation_class": reference.methylation_class,
            "description": reference.description,
            "id": reference.specimens,
            "x": umap_2d[:, 0],
            "y": umap_2d[:, 1],
        })

    logger.info("UMAP done. Matrix created.")
    return (methyl_mtx, umap_df)

def pie_chart(umap_data):
    """Returns plotly pie chart of the methylation classes of the nearest UMAP
    neighbors (according to euclidean 2d-distance) of the sample.
    """
    umap_neighbors = umap_data.umap_df.sort_values(by="distance")[
        1 : UMAP_PLOT_TOP_MATCHES + 1
    ]

    num_per_class = (
        umap_neighbors.groupby(["methylation_class"])
        .size()
        .reset_index(name="counts")
    )
    sample = umap_data.sample
    reference = umap_data.reference
    plot = px.pie(
        num_per_class,
        values="counts",
        names="methylation_class",
        color="methylation_class",
        color_discrete_map=discrete_colors(num_per_class.methylation_class),
        title=(
            f"Nearest UMAP neighbors for {umap_data.sample.name} <br><sup>"
            f"Reference: {reference.name} "
            f"({len(reference.specimens)}"
            f"cases), {len(sample.cpg_overlap)} CpGs</sup>"
        ),
        template="simple_white",
    )
    return plot

class UMAPData:
    """UMAP data container and methods for invoking UMAP plot algorithm."""
    def __init__(self, reference, sample):
        self.sample = sample
        self.reference = reference
        self.methyl_overlap = None
        self.umap_df = None
        self.cu_umap_df = None
        self.plot = None
        self.plot_json = None
        self.cu_plot = None
        self.cu_plot_json = None
        self.pie_chart = None

    @classmethod
    def from_names(cls, reference_name, sample_name):
        return cls(Reference(reference_name), Sample(sample_name))

    @classmethod
    def from_cpgs(cls, reference_name, cpgs):
        sample = Sample.from_cpgs(cpgs)
        return cls(Reference(reference_name), sample)

    def path(self, ending):
        """Returns generic path with corresponding ending."""
        return composite_path(
            NANODIP_REPORTS,
            self.sample.name, self.reference.name, ENDING[ending],
        )

    def make_umap_plot(self):
        """Invoke UMAP plot algorithm and save files to disk."""
        self.methyl_overlap, self.umap_df = dimension_reduction(
            self.reference, self.sample
        )
        self.draw_scatter_plot()
        # If true, sample methylation is part of analysis.
        add_sample = not self.sample.cpgs_only()
        if add_sample:
            self.draw_cu_scatter_plot()
            self.draw_pie_chart()
        self.save_to_disk("all" if add_sample else "no_sample")

    def draw_pie_chart(self):
        """Draw pie chart of nearest UMAP neighbors."""
        self.pie_chart = pie_chart(self)

    def draw_scatter_plot(self):
        """Draws UMAP scatter plot with close-up plot from data."""
        self.plot = umap_plot_from_data(
            self.sample,
            self.reference,
            self.umap_df,
            close_up=False,
        )
        logger.info("UMAP plot generated.")
        # Convert to json.
        self.plot_json = self.plot.to_json()

    def draw_cu_scatter_plot(self):
        """Draws UMAP scatter close-up plot from data."""
        self.cu_umap_df = self.umap_df.sort_values(
            by="distance"
        )[:UMAP_PLOT_TOP_MATCHES + 1]
        self.cu_plot = umap_plot_from_data(
            self.sample,
            self.reference,
            self.cu_umap_df,
            close_up=True,
        )
        logger.info("UMAP close-up plot generated.")
        # Convert to json.
        self.cu_plot_json = self.cu_plot.to_json()

    def save_to_disk(self, params="all"):
        """Saves relevant data to disk."""
        if params == "all":
            obj = ["mmtx", "umtx", "plt", "cplt", "pie", "rnk"]
        elif params == "no_sample":
            obj = ["mmtx", "umtx", "plt"]
        # Save methylation matrix.
        if "mmtx" in obj:
            np.save(self.path("methoverl_npy"), self.methyl_overlap)

        # Save UMAP Matrix.
        if "umtx" in obj:
            self.umap_df.to_csv(self.path("umap_csv"), index=False)

        # Write UMAP plot to disk.
        if "plt" in obj:
            self.plot.write_html(
                self.path("umapall_html"), config=dict({"scrollZoom": True})
            )
            self.plot.write_json(self.path("umapall_json"))
            self.plot.write_image(self.path("umapall_png")) # Time consumption 1.8s

        # Write UMAP close-up plot to disk.
        if "cplt" in obj:
            self.cu_plot.write_html(
                self.path("umaptop_html"), config=dict({"scrollZoom": True})
            )
            self.cu_plot.write_json(self.path("umaptop_json"))
            self.cu_plot.write_image(
                self.path("umaptop_png"), width=600, scale=3
            ) # Time consumption 0.9s

        # Write pie chart to disk.
        if "pie" in obj:
            self.pie_chart.write_image(self.path("pie_png"), width=600, scale=3)

        # Save close up ranking report.
        if "rnk" in obj:
            self.save_ranking_report() # Time consumption 0.4s

    def save_ranking_report(self):
        """Save pdf containing the nearest neighbours from umap analyis."""
        rows = [row for _, row in self.cu_umap_df.iterrows()]
        html_report = render_template("umap_report.html", rows=rows)
        convert_html_to_pdf(html_report, self.path("ranking_pdf"))
        with open(self.path("cpg_cnt"), "w") as f:
            f.write(f"{len(self.sample.cpg_overlap)}")

    def files_on_disk(self, close_up=True):
        """Check if files are on disk."""
        return (
            os.path.exists(self.path("methoverl_npy")) and
            os.path.exists(self.path("umapall_json")) and
            (os.path.exists(self.path("umaptop_json")) or not close_up) and
            os.path.exists(self.path("umap_csv"))
        )

    def read_from_disk(self):
        """Read plot data from disk."""
        # Read UMAP plot as json.
        with open(self.path("umapall_json"), "r") as f:
            self.plot_json = f.read()

        # Read UMAP close-up plot as json if it exists.
        if os.path.exists(self.path("umaptop_json")):
            with open(self.path("umaptop_json"), "r") as f:
                self.cu_plot_json = f.read()

        # Read Methylation Matrix.
        self.methyl_overlap = np.load(self.path("methoverl_npy"), allow_pickle=True)

        # Read UMAP Matrix.
        self.umap_df = pd.read_csv(self.path("umap_csv"))

    def read_precalculated_umap_matrix(self, umap_matrix):
        """Reads precalculated UMAP matrix from disk."""
        path_xlsx = composite_path(
            NANODIP_REPORTS,
            self.sample.name,
            umap_matrix.replace(".xlsx", ""),
            ENDING["umap_xlsx"],
        )
        precalculated_umap = pd.read_excel(
            path_xlsx,
            header=0,
            names=["id", "x", "y"],
            engine="openpyxl",
        )   # TODO better use csv. Time consumption 4.4s
        reference_df = pd.DataFrame(
            zip(
                self.reference.specimens,
                self.reference.methylation_class,
                self.reference.description,
            ),
            columns=["id", "methylation_class", "description"],
        )
        if not self.sample.name in reference_df.id.values:
            reference_df.loc[len(reference_df.index)] = [
                self.sample.name, self.sample.name, "Analysis Sample"
            ]
        self.umap_df = pd.merge(precalculated_umap, reference_df, on = "id")
        umap_sample = self.umap_df[["x", "y"]].loc[
            self.umap_df.id==self.sample.name
        ].values
        self.umap_df["distance"] = [
            np.linalg.norm([z.x, z.y] - umap_sample)
            for z in self.umap_df.itertuples()
        ]
        self.umap_df = self.umap_df.sort_values(by="distance")
