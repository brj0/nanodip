"""
## Plots

Functions for creating methylation UMAP plot.
Functions for creating Copy Number Variation plot.
"""

# start_external_modules
from plotly.io import write_json, from_json
import bisect
import csv
import logging
import math
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import time
# end_external_modules

# start_internal_modules
from config import (
    CNV_GRID,
    CNV_LINK,
    ENDING,
    NANODIP_REPORTS,
    PLOTLY_RENDER_MODE,
    UMAP_PLOT_TOP_MATCHES,
)
from data import (
    Reference,
    Genome,
    Sample,
    get_reference_methylation,
    get_sample_methylation,
)
from utils import (
    convert_html_to_pdf,
    discrete_colors,
    render_template,
    composite_path,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)

def umap_plot_from_data(sample, reference, umap_data_frame, close_up):
    """Create and return umap plot from UMAP data.

    Args:
        sample: sample data
        reference: reference data
        umap_data_frame: pandas data frame containing umap info. First
            row corresponds to sample.
        close_up: bool to indicate if only top matches should be plotted.
    """
    umap_sample = umap_data_frame.iloc[0]
    umap_title = (
        f"UMAP for {sample.name} <br><sup>Reference: {reference.name} "
        f"({len(reference.specimens)} cases), "
        f"{len(sample.cpg_overlap)} CpGs </sup>"
    )
    if close_up:
        umap_title = "Close-up " + umap_title
    methyl_classes = umap_data_frame.methylation_class[1:].to_list()
    methyl_classes.sort()
    umap_plot = px.scatter(
        umap_data_frame,
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

    # If close-up add hyperlinks for all references and draw circle
    if close_up:
        umap_plot.update_traces(marker=dict(size=5))
        # Add hyperlinks
        for _, row in umap_data_frame.iloc[1:].iterrows():
            URL = CNV_LINK % row["id"]
            umap_plot.add_annotation(
                x=row["x"],
                y=row["y"],
                text=f"<a href='{URL}' target='_blank'>&nbsp;</a>",
                showarrow=False,
                arrowhead=1,
            )
        # Draw circle
        radius = umap_data_frame["distance"].iloc[-1]
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

def umap_data_frame(sample, reference):
    """Create UMAP methylation analysis matrix.

    Args:
        sample: sample to analyse
        reference: reference data
    """
    import umap #Moved here due to long loading time (13.5s)

    logger.info(f"Start UMAP for {sample.name} / {reference.name}.")
    logger.info(reference)

    if sample.methyl_df.empty:
        logger.info("UMAP done. No Matrix created, no overlapping data.")
        raise ValueError("Sample has no overlapping CpG's with reference.")

    # Calculate overlap of sample CpG's with reference CpG's (some probes have
    # been skipped from the reference set, e.g. sex chromosomes).
    sample.set_cpg_overlap(reference)
    logger.info(sample)

    # Extract reference and sample methylation according to CpG overlap.
    reference_methylation = get_reference_methylation(sample,
                                                      reference)
    logger.info(f"Reference methylation extracted:\n{reference_methylation}")

    sample_methylation = get_sample_methylation(sample, reference)
    logger.info(f"Sample methylation extracted:\n{sample_methylation}")

    # Calculate UMAP Nx2 Matrix. Time intensive (~1min).
    methyl_overlap = np.vstack([sample_methylation, reference_methylation])
    logger.info("UMAP algorithm initiated.")
    umap_2d = umap.UMAP(verbose=True).fit_transform(methyl_overlap)
    logger.info("UMAP algorithm done.")

    # Free memory
    del reference_methylation
    del sample_methylation

    umap_sample = umap_2d[0]
    umap_df = pd.DataFrame({
        "distance": [np.linalg.norm(z - umap_sample) for z in umap_2d],
        "methylation_class": [sample.name] + reference.methylation_class,
        "description":  ["Analysis sample"] + reference.description,
        "id": [sample.name] + reference.specimens,
        "x": umap_2d[:,0],
        "y": umap_2d[:,1],
    })

    logger.info("UMAP done. Matrix created.")

    return (methyl_overlap, umap_df)

def pie_chart(umap_data):
    """Create pie chart as png of the nearst umap neighbors (according to 
    euclidean 2d-distance).
    """
    umap_neighbors = umap_data.umap_df.sort_values(by="distance")[
        1 : UMAP_PLOT_TOP_MATCHES + 1
    ]

    num_per_class = (
        umap_neighbors.groupby(["methylation_class"])
        .size()
        .reset_index(name="counts")
    )

    if sum(num_per_class.counts) < UMAP_PLOT_TOP_MATCHES:
        raise ValueError("Not enough reference values.")
    sample = umap_data.sample
    reference = umap_data.reference

    pie_chart = px.pie(
        num_per_class,
        values="counts",
        names="methylation_class",
        color="methylation_class",
        color_discrete_map=discrete_colors(num_per_class.methylation_class),
        title=(
            f"Nearest UMAP neighbors for {umap_data.sample_name} <br><sup>"
            f"Reference: {reference.name} "
            f"({len(reference.specimens)}"
            f"cases), {len(sample.cpg_overlap)} CpGs</sup>"
        ),
        template="simple_white",
    )
    return pie_chart

def get_bin_edges(n_bins, genome):
    """Returns sequence of {n_bin} equal sized bins on chromosomes. Every bin is
    limited to one chromosome."""
    edges = np.linspace(0, len(genome), num=n_bins + 1).astype(int)
    # limit bins to only one chromosome
    for chrom_edge in genome.chrom.offset:
        i_nearest = np.abs(edges - chrom_edge).argmin()
        edges[i_nearest] = chrom_edge
    return edges

def get_cnv(read_positions, genome):
    """Return CNV."""
    expected_reads_per_bin = 30
    n_bins = len(read_positions)//expected_reads_per_bin
    read_start_positions = [i[0] for i in read_positions]

    copy_numbers, bin_edges = np.histogram(
        read_start_positions,
        bins=get_bin_edges(n_bins, genome),
        range=[0, len(genome)],
    )

    bin_midpoints = (bin_edges[1:] + bin_edges[:-1])/2
    expected_reads = np.diff(bin_edges) * len(read_positions) / len(genome)
    cnv = [(x - e)/e for x, e in zip(copy_numbers, expected_reads)]
    return bin_midpoints, copy_numbers

def cnv_grid(genome):
    """Makes chromosome grid layout for CNV Plot and saves it on disk. If
    available grid is read from disk.
    """
    # Check if grid exists and return if available.
    grid_path = CNV_GRID
    if os.path.exists(grid_path):
        with open(grid_path, "r") as f:
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
        grid.add_vline(x=i, line_color="black",
                           line_dash="dot", line_width=1)
    # Vertical line: shromosomes.
    for i in genome.chrom.offset.tolist() + [len(genome)]:
        grid.add_vline(x=i, line_color="black", line_width=1)
    # Save to disk
    grid.write_json(grid_path)
    return grid

def cnv_plot_from_data(data_x, data_y, E_y, sample_name, read_num, genome):
    """Create CNV plot from CNV data.

    Args:
        data_x: x-Values to plot.
        data_y: y-Values to plot.
        E_y: expected y-Value.
        sample_name: Name of sample.
        read_num: Number of read reads.
        genome: Reference Genome.
    """
    grid = cnv_grid(genome)
    # Expected value: horizontal line.
    grid.add_hline(y=E_y, line_color="black", line_width=1)
    plot = px.scatter(
        x=data_x,
        y=data_y,
        labels={
            "x":f"Number of mapped reads: {read_num}",
            "y":f"Copy numbers per {round(len(genome)/(len(data_x) * 1e6), 2)} MB"
        },
        title=f"Sample ID: {sample_name}",
        color=data_y,
        range_color=[E_y*0, E_y*2],
        color_continuous_scale="Portland",
        render_mode=PLOTLY_RENDER_MODE,
    )
    plot.update_traces(
        hovertemplate="Copy Numbers = %{y} <br>",
    )
    plot.update_layout(
        grid.layout,
        yaxis_range = [-0.5, 2*E_y],
    )
    return plot

def number_of_reads(read_start_pos, interval):
    """Return the number of starting sequences whithin interval. Reads must
    be sorted in ascending order."""
    left, right = interval
    i_left = bisect.bisect_left(read_start_pos, left)
    i_right = bisect.bisect_left(read_start_pos, right)
    return len(read_start_pos[i_left:i_right])

def cnv_plot(sample, bin_midpoints, cnv, genome):
    """Create a genome-wide copy number plot and save data on dist."""
    logger.info(f"CNVP start")
    logger.info(f"Read positions:\n{sample.reads[:100]}")
    logger.info(f"Bin midpoints:\n{bin_midpoints}")
    logger.info(f"CNV:\n{cnv}")

    avg_read_per_bin = len(sample.reads) // len(bin_midpoints)

    plot = cnv_plot_from_data(
        data_x=bin_midpoints,
        data_y=cnv,
        E_y = avg_read_per_bin,
        sample_name=sample.name,
        read_num=len(sample.reads),
        genome=genome,
    )

    logger.info(f"CNVP done")
    return plot

class UMAPData:
    """Umap data container and methods for invoking umap plot algorithm."""
    def __init__(self, sample_name, reference_name):
        self.sample_name = sample_name
        self.reference_name = reference_name
        self.sample = Sample(self.sample_name)
        self.reference = Reference(self.reference_name)

    def make_umap_plot(self):
        """Invoke umap plot algorithm and save to disk."""
        self.methyl_overlap, self.umap_df = umap_data_frame(
            self.sample, self.reference
        )
        # Save Methylation Matrix.
        np.save(self.path("methyl"), self.methyl_overlap)

        self.draw_scatter_plots()
        self.draw_pie_chart()

        # Save close up ranking report.
        self.save_ranking_report() # Time consumption 0.4s

    def draw_scatter_plots(self):
        self.plot = umap_plot_from_data(
            self.sample,
            self.reference,
            self.umap_df,
            close_up=False,
        )
        logger.info("UMAP plot generated.")
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
        self.plot_json = self.plot.to_json()
        self.cu_plot_json = self.cu_plot.to_json()

        # Save UMAP Matrix.
        self.umap_df.to_csv(self.path("umap_csv"), index=False)

        # Write UMAP plot to disk.
        file_path = self.path("umap_all")
        self.plot.write_html(file_path, config=dict({"scrollZoom": True}))
        self.plot.write_json(file_path[:-4] + "json")
        self.plot.write_image(file_path[:-4] + "png") # Time consumption 1.8s

        # Write UMAP close-up plot to disk.
        file_path = self.path("umap_top")
        self.cu_plot.write_html(file_path, config=dict({"scrollZoom": True}))
        self.cu_plot.write_json(file_path[:-4] + "json")

        # Time consumption 0.9s
        self.cu_plot.write_image(
            file_path[:-4] + "png", width=450, height=400, scale=3
        )

    def draw_pie_chart(self):
        self.pie_chart = pie_chart(self)

        # Write pie chart to disk.
        self.pie_chart.write_image(
            self.path("pie"), width=450, height=400, scale=3
        )

    def save_ranking_report(self):
        """Save pdf containing the nearest neighbours from umap analyis."""
        rows = [row for _, row in self.cu_umap_df.iterrows()]

        html_report = render_template("umap_report.html", rows=rows)

        convert_html_to_pdf(html_report, self.path("ranking"))

        with open(self.path("cpg_cnt"), "w") as f:
            f.write("%s" % len(self.sample.cpg_overlap))

    def path(self, ending):
        return composite_path(
            NANODIP_REPORTS,
            self.sample_name, self.reference_name, ENDING[ending],
        )

    def files_on_disk(self):
        return (
            os.path.exists(self.path("methyl")) and
            os.path.exists(self.path("umap_all_json")) and
            os.path.exists(self.path("umap_top_json")) and
            os.path.exists(self.path("umap_csv"))
        )

    def read_from_disk(self):
        methyl_overlap_path = self.path("methyl")
        plot_path = self.path("umap_all_json")
        cu_plot_path = self.path("umap_top_json")

        # Read UMAP plot as json.
        with open(plot_path, "r") as f:
            self.plot_json = f.read()

        # Read UMAP close-up plot as json.
        with open(cu_plot_path, "r") as f:
            self.cu_plot_json = f.read()

        # Read Methylation Matrix.
        self.methyl_overlap = np.load(methyl_overlap_path,
            allow_pickle=True)

        # Read UMAP Matrix.
        self.umap_df = pd.read_csv(self.path("umap_csv"))

    def read_precalculated_umap_matrix(self, umap_matrix):
        self.reference = Reference(self.reference_name)
        path_xlsx = composite_path(
            NANODIP_REPORTS,
            self.sample_name,
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
        if not self.sample_name in reference_df.id.values:
            reference_df.loc[len(reference_df.index)] = [
                self.sample_name, self.sample_name, "Analysis Sample"
            ]
        self.umap_df = pd.merge(
            precalculated_umap, reference_df, on = "id"
        )
        umap_sample = self.umap_df[["x", "y"]].loc[
            self.umap_df.id==self.sample_name
        ].values
        self.umap_df["distance"] = [
            np.linalg.norm([z.x, z.y] - umap_sample)
            for z in self.umap_df.itertuples()
        ]
        self.umap_df = self.umap_df.sort_values(by="distance")

class CNVData:
    """CNV data container and methods for invoking cnv plot algorithm."""
    genome = Genome()
    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.base_path = composite_path(NANODIP_REPORTS, sample_name, "")

    def read_from_disk(self):
        plot_path = self.base_path + ENDING["cnv_json"]
        genes_path = self.base_path + ENDING["genes"]
        with open(plot_path, "r") as f:
            self.plot_json = f.read()
        self.plot = from_json(self.plot_json)
        self.genes = pd.read_csv(genes_path)

    def files_on_disk(self):
        plot_path = self.base_path + ENDING["cnv_json"]
        genes_path = self.base_path + ENDING["genes"]
        return (
            os.path.exists(plot_path) and
            os.path.exists(genes_path)
        )

    def make_cnv_plot(self):
        self.sample = Sample(self.sample_name)
        self.sample.set_reads() # time consumption 2.5s
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
        self.genes = self.gene_cnv(
            len(CNVData.genome) // len(self.bin_midpoints)
        )
        self.relevant_genes = self.genes.loc[self.genes.relevant]
        self.save_to_disk()

    def save_to_disk(self):
        self.plot.write_html(
            self.base_path + ENDING["cnv_html"],
            config=dict({"scrollZoom": True}),
        )
        write_json(self.plot, self.base_path + ENDING["cnv_json"])
        # time consuming operation (1.96s)
        self.plot.write_image(
            self.base_path + ENDING["cnv_png"], width=1280, height=720, scale=3,
        )
        with open(self.base_path + ENDING["aligned_reads"], "w") as f:
            f.write("%s" % len(self.sample.reads))
        with open(self.base_path + ENDING["reads_csv"], "w") as f:
            write = csv.writer(f)
            write.writerows(self.sample.reads)
        self.genes.to_csv(self.base_path + ENDING["genes"], index=False)
        self.relevant_genes.to_csv(
            self.base_path + ENDING["relevant_genes"],
            index=False,
        )

    def gene_cnv(self, bin_width):
        genes = CNVData.genome.genes
        genes["interval"] = list(zip(genes.start, genes.end))
        read_start_pos = [i[0] for i in self.sample.reads]
        read_start_pos.sort()
        num_reads = len(read_start_pos)
        genes["cn_obs"] = genes.interval.apply(
            lambda z: number_of_reads(read_start_pos, z)
        )
        genes["cn_per_mega_base"] = genes.apply(
            lambda z: z["cn_obs"]/z["len"] * 1e6, # TODO auto draw extreme values
            axis=1,
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
        gene_pos = self.genes.loc[self.genes.name.isin(genes)]
        return gene_pos

    def plot_cnv_and_genes(self, gene_names):
        genes = self.get_gene_positions(gene_names)
        plot = go.Figure(self.plot)
        plot.add_trace(
            go.Scatter(
                customdata=genes[[
                    "name",          # 0
                    "loc",           # 1
                    "transcript",    # 2
                    "len",           # 3
                ]],
                hovertemplate=(
                    "Copy numbers = %{y} <br>"
                    "<b> %{customdata[0]} </b> <br>"
                    "%{customdata[1]} "
                    "(hg19 %{customdata[2]}) <br>"
                    "%{customdata[3]} bases <br>"
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
                y=genes.cn_obs,
            ))
        return plot.to_json()
