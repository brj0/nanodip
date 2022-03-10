"""
## Plots

Create Methylation UMAP plot.
Create Copy Number Variation plot.
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
    CNV_URL_PREFIX,
    CNV_URL_SUFFIX,
    ENDINGS,
    NANODIP_REPORTS,
    PLOTLY_RENDER_MODE,
    UMAP_PLOT_TOP_MATCHES,
)
from data import (
    ReferenceData,
    ReferenceGenome,
    SampleData,
    get_reference_methylation,
    get_sample_methylation,
)
from utils import (
    convert_html_to_pdf,
    render_template,
)
# end_internal_modules


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
    umap_title = f"UMAP for {sample.name} against {reference.name}, "\
        + f"{len(reference.annotated_specimens)} reference cases, "\
        + f"{len(sample.cpg_overlap)} CpGs"
    if close_up:
        umap_title = "Close-up " + umap_title
    umap_plot = px.scatter(
        umap_data_frame,
        x="x",
        y="y",
        labels={"x":"UMAP 0", "y":"UMAP 1", "color":"WHO class"},
        title=umap_title,
        color="methylation_class",
        hover_name="id",
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
            umap_plot.add_annotation(
                x=row["x"],
                y=row["y"],
                text="<a href='" + CNV_URL_PREFIX + row["id"]
                    + CNV_URL_SUFFIX
                    + "' target='_blank'>&nbsp;</a>",
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
    import umap #TODO move to beginning

    logger.info(
        f"UMAP Plot initiated for {sample.name} and reference {reference.name}."
    )
    logger.info(f"Reference Annotation:\n{reference.annotation}")
    logger.info(f"Reference CpG Sites No:\n{len(reference.cpg_sites)}")
    logger.info(f"Reference Specimens No:\n{len(reference.specimens)}")
    logger.info(
        f"Reference Annotated specimens: {len(reference.annotated_specimens)}"
    )
    logger.info(f"Sample CpG Sites No:\n{len(sample.cpg_sites)}")
    logger.info(f"Sample CpG overlap No before:\n{sample.cpg_overlap}")

    # Calculate overlap of sample CpG's with reference CpG's (some probes have
    # been skipped from the reference set, e.g. sex chromosomes).
    sample.set_cpg_overlap(reference)
    logger.info(f"Sample read. CpG overlap No after:\n{len(sample.cpg_overlap)}")

    if not sample.cpg_overlap:
        logger.info("UMAP done. No Matrix created, no overlapping data.")
        raise ValueError("Sample has no overlapping CpG's with reference.")

    # Extract reference and sample methylation according to CpG overlap.
    reference_methylation = get_reference_methylation(sample,
                                                      reference)
    logger.info(f"""Reference methylation extracted:
                {reference_methylation}""")
    sample_methylation = get_sample_methylation(sample, reference)
    logger.info(f"""Sample methylation extracted:
                {sample_methylation}""")
    logger.info("UMAP algorithm initiated.")

    # Calculate UMAP Nx2 Matrix. Time intensive (~1min).
    methyl_overlap = np.vstack([sample_methylation, reference_methylation])
    umap_2d = umap.UMAP(verbose=True).fit_transform(methyl_overlap)

    # Free memory
    del reference_methylation
    del sample_methylation

    logger.info("UMAP algorithm done.")

    umap_sample = umap_2d[0]
    umap_df = pd.DataFrame({
        "distance": [np.linalg.norm(z - umap_sample) for z in umap_2d],
        "methylation_class":  [sample.name] + reference.methylation_class,
        "description":  ["Sample tested"] + reference.description,
        "id": [sample.name] + reference.specimen_ids,
        "x": umap_2d[:,0],
        "y": umap_2d[:,1],
    })

    logger.info("UMAP done. Matrix created.")

    return (methyl_overlap, umap_df)

def get_bin_edges(n_bins, genome):
    """Returns sequence of {n_bin} equal sized bins on chromosomes. Every bin is
    limited to one chromosome."""
    if n_bins < 100:
        raise ValueError("Binwidth too small.")
    edges = np.linspace(0, genome.length, num=n_bins + 1).astype(int)
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
        read_start_positions, bins=get_bin_edges(n_bins, genome), range=[0, genome.length]
    )

    bin_midpoints = (bin_edges[1:] + bin_edges[:-1])/2
    expected_reads = np.diff(bin_edges) * len(read_positions) / genome.length
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
            range=[0, genome.length],
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
        width=1700,
        height=900,
        template="simple_white",
    )
    # Vertical line: centromere.
    for i in genome.chrom.centromere_offset:
        grid.add_vline(x=i, line_color="black",
                           line_dash="dot", line_width=1)
    # Vertical line: shromosomes.
    for i in genome.chrom.offset.tolist() + [genome.length]:
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
    cnv_plot = px.scatter(
        x=data_x,
        y=data_y,
        labels={
            "x":f"Number of mapped reads: {read_num}",
            "y":f"Copy numbers per {round(genome.length/(len(data_x) * 1e6), 2)} MB"
        },
        title=f"Sample ID: {sample_name}",
        color=data_y,
        range_color=[E_y*0, E_y*2],
        color_continuous_scale="Portland",
        render_mode=PLOTLY_RENDER_MODE,
    )
    cnv_plot.update_traces(
        hovertemplate="Copy Numbers = %{y} <br>",
    )
    cnv_plot.update_layout(
        grid.layout,
        yaxis_range = [-0.5, 2*E_y],
    )
    return cnv_plot

def number_of_reads(read_start_pos, interval):
    """Return the number of starting sequences whithin interval. Reads must
    be sorted in ascending order."""
    left, right = interval
    i_left = bisect.bisect_left(read_start_pos, left)
    i_right = bisect.bisect_left(read_start_pos, right)
    return len(read_start_pos[i_left:i_right])

def get_cnv_plot(sample, bin_midpoints, cnv, genome):
    """Create a genome-wide copy number plot and save data on dist."""
    logger.info(f"CNVP start")
    logger.info(f"Read positions:\n{sample.reads[:100]}")

    logger.info(f"Bin midpoints:\n{bin_midpoints}")
    logger.info(f"CNV:\n{cnv}")

    avg_read_per_bin = len(sample.reads) // len(bin_midpoints)

    cnv_plot = cnv_plot_from_data(
        data_x=bin_midpoints,
        data_y=cnv,
        E_y = avg_read_per_bin,
        sample_name=sample.name,
        read_num=len(sample.reads),
        genome=genome,
    )

    logger.info(f"CNVP done")
    return cnv_plot

class UMAPData:
    """Umap data container and methods for invoking umap plot algorithm."""
    def __init__(self, sample_name, reference_name):
        self.sample_name = sample_name
        self.reference_name = reference_name
        self.path = os.path.join(
            NANODIP_REPORTS,
            f"{sample_name}_{reference_name}",
        )

    def make_umap_plot(self):
        """Invoke umap plot algorithm and save to disk."""
        self.sample = SampleData(self.sample_name)
        self.reference = ReferenceData(self.reference_name) # time: 3.6s
        self.methyl_overlap, self.umap_df = umap_data_frame(
            self.sample, self.reference
        )
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

        self.save_to_disk()

    def save_ranking_report(self):
        """Save pdf containing the nearest neighbours from umap analyis."""
        rows = [row for _, row in self.cu_umap_df.iterrows()]

        html_report = render_template("umap_report.html", rows=rows)

        file_path = os.path.join(
            NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         ENDINGS["ranking"],
                        ),
        )
        convert_html_to_pdf(html_report, file_path)

        file_path = os.path.join(
            NANODIP_REPORTS,
            "%s%s" % (self.sample_name, ENDINGS["cpg_cnt"]),
        )
        with open(file_path, "w") as f:
            f.write("%s" % len(self.sample.cpg_overlap))

    def save_to_disk(self):
        # Save Methylation Matrix.
        file_path = os.path.join(NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         ENDINGS["methyl"])
        )
        np.save(file_path, self.methyl_overlap)

        # Save UMAP Matrix.
        file_path = os.path.join(NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         ENDINGS["umap_csv"])
        )
        self.umap_df.to_csv(file_path, index=False)

        # Write UMAP plot to disk.
        file_path = os.path.join(
            NANODIP_REPORTS,
            f"%s_%s%s" % (self.sample_name, self.reference_name,
                          ENDINGS["umap_all"]),
        )
        self.plot.write_html(file_path, config=dict({"scrollZoom": True}))
        self.plot.write_json(file_path[:-4] + "json")
        self.plot.write_image(file_path[:-4] + "png") # Time consumption 1.8s

        # Write UMAP close-up plot to disk.
        file_path = os.path.join(
            NANODIP_REPORTS,
            f"%s_%s%s" % (self.sample_name, self.reference_name,
                          ENDINGS["umap_top"]),
        )
        self.cu_plot.write_html(file_path, config=dict({"scrollZoom": True}))
        self.cu_plot.write_json(file_path[:-4] + "json")
        self.cu_plot.write_image(file_path[:-4] + "png") # Time consumption 0.9s

        # Save close up ranking report.
        self.save_ranking_report() # Time consumption 0.4s


    def files_on_disk(self):
        methyl_overlap_path = os.path.join(NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         ENDINGS["methyl"])
        )
        plot_path = self.path + ENDINGS["umap_all_json"]
        cu_plot_path = self.path + ENDINGS["umap_top_json"]

        return (os.path.exists(plot_path) and
                os.path.exists(cu_plot_path) and
                os.path.exists(methyl_overlap_path))

    def read_from_disk(self):
        methyl_overlap_path = os.path.join(NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         ENDINGS["methyl"])
        )
        plot_path = self.path + ENDINGS["umap_all_json"]
        cu_plot_path = self.path + ENDINGS["umap_top_json"]

        # Read UMAP plot as json.
        with open(plot_path, "r") as f:
            self.plot_json = f.read()
        #self.plot = from_json(self.plot_json)

        # Read UMAP close-up plot as json.
        with open(cu_plot_path, "r") as f:
            self.cu_plot_json = f.read()
        #self.cu_plot = from_json(self.cu_plot_json)

        # Read Methylation Matrix.
        self.methyl_overlap = np.load(methyl_overlap_path,
            allow_pickle=True)


def TODO():
    path = umap_output_path(sample, reference, close_up=True)
    with open(path["html"], "w") as f:
        f.write("<html><body>No data to plot.</body></html>")


class CNVData:
    """CNV data container and methods for invoking cnv plot algorithm."""
    genome = ReferenceGenome()
    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.path = os.path.join(NANODIP_REPORTS, f"{sample_name}")

    def read_from_disk(self):
        plot_path = self.path + ENDINGS["cnv_json"]
        genes_path = self.path + ENDINGS["genes"]
        with open(plot_path, "r") as f:
            self.plot_json = f.read()
        self.plot = from_json(self.plot_json)
        self.genes = pd.read_csv(genes_path)

    def files_on_disk(self):
        plot_path = self.path + ENDINGS["cnv_json"]
        genes_path = self.path + ENDINGS["genes"]
        return (
            os.path.exists(plot_path) and
            os.path.exists(genes_path)
        )

    def make_cnv_plot(self):
        self.sample = SampleData(self.sample_name)
        self.sample.set_reads() # time consumption 2.5s
        self.bin_midpoints, self.cnv = get_cnv(
            self.sample.reads,
            CNVData.genome,
        )
        self.plot = get_cnv_plot(
            sample=self.sample,
            bin_midpoints=self.bin_midpoints,
            cnv=self.cnv,
            genome=CNVData.genome,
        )
        self.plot_json = self.plot.to_json()
        self.genes = self.gene_cnv(
            CNVData.genome.length // len(self.bin_midpoints)
        )
        self.relevant_genes = self.genes.loc[self.genes.relevant]
        self.save_to_disk()

    def save_to_disk(self):
        self.plot.write_html(
            self.path + ENDINGS["cnv_html"],
            config=dict({"scrollZoom": True}),
        )
        write_json(self.plot, self.path + ENDINGS["cnv_json"])
        if not os.path.exists(self.path + ENDINGS["cnv_png"]):
            # time consuming operation (1.96s)
            self.plot.write_image(
                self.path + ENDINGS["cnv_png"], width=1280, height=720
            )
        with open(self.path + ENDINGS["aligned_reads"], "w") as f:
            f.write("%s" % len(self.sample.reads))
        with open(self.path + ENDINGS["reads_csv"], "w") as f:
            write = csv.writer(f)
            write.writerows(self.sample.reads)
        self.genes.to_csv(self.path + ENDINGS["genes"], index=False)
        self.relevant_genes.to_csv(
            self.path + ENDINGS["relevant_genes"],
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
        genes["cn_norm"] = genes.apply(
            lambda z: z["cn_obs"]/z["len"] * bin_width, # TODO exp/norm not compatible
            axis=1,
        )
        genes["cn_exp"] = genes.apply(
            lambda z: len(self.sample.reads)*z["len"]/CNVData.genome.length,
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


# create and cofigure logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(levelname)s %(asctime)s %(lineno)d - %(message)s")
file_handler = logging.FileHandler(
    os.path.join(NANODIP_REPORTS, "nanodip.log"),
    "w",
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

## output to console
# stream_handler = logging.StreamHandler()
# logger.addHandler(stream_handler)
