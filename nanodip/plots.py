from plotly.io import write_json, from_json
from utils import render_template
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
import xhtml2pdf.pisa

import config
import data


def convert_html_to_pdf(source_html, output_file):
    """Create PDF from html-string."""
    with open(output_file, "w+b") as f:
        pisa_status = xhtml2pdf.pisa.CreatePDF(source_html, dest=f)
    return pisa_status.err

def umap_output_path(sample, reference, close_up):
    """Generate path to save umap plots in html and png."""
    ending = "top" if close_up else "all"
    file_path = {}
    for i in ["html", "png"]:
        file_path[i] = os.path.join(
            config.NANODIP_REPORTS,
            f"{sample.name}_{reference.name}_UMAP_{ending}.{i}",
        )
    return file_path

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
        x=umap_data_frame["x"],
        y=umap_data_frame["y"],
        labels={"x":"UMAP 0", "y":"UMAP 1", "color":"WHO class"},
        title=umap_title,
        color=umap_data_frame["methylation_class"],
        hover_name=umap_data_frame["id"],
        render_mode=config.PLOTLY_RENDER_MODE,
    )
    umap_plot.add_annotation(
        x=umap_sample["x"],
        y=umap_sample["y"],
        text=sample.name,
        showarrow=True,
        arrowhead=1,
    )
    umap_plot.update_yaxes(scaleanchor = "x", scaleratio = 1)
    # If close-up add hyperlinks for all references and draw circle
    if close_up:
        umap_plot.update_traces(marker=dict(size=5))
        # Add hyperlinks
        for _, row in umap_data_frame.iloc[1:].iterrows():
            umap_plot.add_annotation(
                x=row["x"],
                y=row["y"],
                text="<a href='" + config.CNV_URL_PREFIX + row["id"]
                    + config.CNV_URL_SUFFIX
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
            line_color="black",
            line_width=0.5,
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
    reference_methylation = data.get_reference_methylation(sample,
                                                           reference)
    logger.info(f"""Reference methylation extracted:
                {reference_methylation}""")
    sample_methylation = data.get_sample_methylation(sample, reference)
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
        "id": [sample.name] + reference.specimen_ids,
        "x": umap_2d[:,0],
        "y": umap_2d[:,1],
    })

    logger.info("UMAP done. Matrix created.")

    return (methyl_overlap, umap_df)

class SortedList:
    """Container of reversely sorted list with efficient pop and min
    lookup function."""
    def __init__(self, list, key=None):
        self.list = list
        self.key = key if key else lambda x: x
        self.list.sort(reverse=True, key=key)
    def __repr__(self):
        return str(self.list)
    def insert(self, value):
        self.list.append(value)
        self.list.sort(reverse=True, key=self.key)
    def extend(self, list):
        for l in list:
            self.insert(l)
    def min(self):
        return self.key(self.list[-1])
    def pop(self):
        """Pop all list values equal to the list min.

        Returns: Poped items.
        """
        list_min = self.min()
        poped = []
        while self.list and self.min() == list_min:
            poped.append(self.list.pop())
        return poped
    def __len__(self):
        # Necessary for boolean tests
        return len(self.list)

def get_copy_numbers(genome_length, read_intervals):
    """Calculate Copy numbers.

    Args:
        genome_length: Length of genome
        read_intervals: List of 2-element lists [a,b] containing left
            read end "a" and read end "b" (more exactly: first element
            after read).
    Returns:
        step_pos: Position of coverage steps
        copy_num: Coverage number.
    """
    remaining = SortedList(read_intervals.copy(), key=lambda x: x[0])
    next_rights = SortedList([])
    copy_num = []
    step_pos = [] if remaining.min() == 0 else [0]
    get_right = lambda x: [y[1] for y in x]
    no_gap_between_reads = False

    while remaining or next_rights:
        left = None if not remaining else remaining.min()
        right = None if not next_rights else next_rights.min()
        # Process the next left read ends.
        if left is not None and (right is None or left < right):
            if no_gap_between_reads:
                no_gap_between_reads = False
            else:
                copy_num.append(len(next_rights))
                step_pos.append(left)
            poped_rights = get_right(remaining.pop())
            next_rights.extend(poped_rights)
        # Process the next right read ends. Set flag for case of
        # seamless alignment.
        if left is not None and right is not None and left == right:
            no_gap_between_reads = True
            copy_num.append(len(next_rights))
            step_pos.append(right)
            next_rights.pop()
        # Process the next right read.
        if right is not None and (left is None or right < left):
            copy_num.append(len(next_rights))
            step_pos.append(right)
            next_rights.pop()
    if step_pos[-1] != genome_length:
        step_pos.append(genome_length)
        copy_num.append(0)
    return step_pos, copy_num

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

def get_cnv_from_read_start(read_positions, genome):
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
    return bin_midpoints, copy_numbers#cnv

def mean_local_coverage(read_positions, genome):
    """Calculate mean coverage whithin evenly spaced bins.

    Args:
        read_intervals: List of read start/end positions.
    Returns:
        bin_midpoints: List of midpoints of bins of genome.
        cnv: Average coverage. Non covered sections will be ignored.
    """
    expected_reads_per_bin = 20
    n_bins = len(read_positions)//expected_reads_per_bin
    step = math.ceil(genome.length / n_bins)
    bin_edges = get_bin_edges(n_bins, genome)
    step_pos, copy_num = get_copy_numbers(genome.length, read_positions)
    cnv = []
    bin_midpoints = []
    total_read_len = sum([i[1] - i[0] for i in read_positions])
    for i in zip(bin_edges[:-1], bin_edges[1:]):
        expected = (i[1] - i[0])/genome.length*total_read_len
        observed = copy_number_sum(step_pos, copy_num, i)
        if observed > 0:
            normalized = (observed - expected) / expected
            cnv.append(normalized)
            bin_midpoints.append(np.mean(i))
    return bin_midpoints, cnv

def cnv_plot_from_data(data_x, data_y, avg_y, sample, read_num, genome):
    """Create CNV plot from CNV data.

    Args:
        data_x: x-Values to plot.
        data_y: y-Values to plot.
        sample: Sample object containing sample name.
        read_num: Number of read reads.
        genome: Reference Genome.
    """
    cnv_plot = px.scatter(
        x=data_x,
        y=data_y,
        labels={"x":f"Number of mapped reads: {read_num}",
            "y":f"Relative CNV per {round(genome.length/(len(data_x) * 1e6), 2)} MB bin"},
        title=f"Sample ID: {sample.name}",
        color=data_y,
        range_color=[0, 2*avg_y],
        template="simple_white",
        #hover_name=bin_midpoints,
        color_continuous_scale="RdBu",
        render_mode=config.PLOTLY_RENDER_MODE,
    )
    cnv_plot.update_layout(coloraxis_showscale=False)
    cnv_plot.update_layout(
        xaxis = dict(
            zeroline=False,
            tickmode="array",
            tickvals=genome.chrom.center,
            ticktext=genome.chrom.name,
            showgrid=False,
            ticks="outside",
            tickson="boundaries",
            ticklen=10
        )
    )
    cnv_plot.add_hline(y=avg_y, line_color="black", line_width=1)
    for i in genome.chrom.centromere_offset:
        cnv_plot.add_vline(x=i, line_color="black",
                           line_dash="dot", line_width=1)
    for i in genome.chrom.offset.tolist() + [genome.length]:
        cnv_plot.add_vline(x=i, line_color="black", line_width=1)
    cnv_plot.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        range=[0, genome.length],
    )
    cnv_plot.update_yaxes(showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        range = [0, 2*avg_y],
        #range=[-1.5, 1.5],
    )
    cnv_plot.update_layout(width=1700, height=900)
    return cnv_plot

def binary_search(search_list, value, upper=False):
    """Return index of biggest value in the search list that is less or equal
    than value, respecively the smallest value greater or equal value.

    Args:
        search_list: ascending ordered list
        value: value to search
        upper: bool to indicate if smallest or biggest value should be
            searched
    Returns:
        Biggest value in search_list <=value if upper and
        smallest value in search_list >=value if not upper.
    """
    if not upper:
        pos = bisect.bisect_left(search_list, value)
        if pos == 0:
            return pos
        elif pos == len(search_list):
            return pos - 1
        elif search_list[pos] == value:
            return pos
        else:
            return pos - 1
    elif upper:
        pos = bisect.bisect_right(search_list, value)
        if pos == 0:
            return pos
        elif pos == len(search_list):
            return pos - 1
        elif search_list[pos - 1] == value:
            return pos - 1
        else:
            return pos

def copy_number_sum(step_pos, copy_numbers, interval):
    """Integrate step function."""
    if not len(interval) == 2:
        raise ValueError("invalid range")
    if not len(step_pos) == len(copy_numbers) + 1:
        raise ValueError("invalid lengths of step_pos/copy_numbers")
    left = binary_search(step_pos, interval[0])
    right = binary_search(step_pos, interval[1], upper=True)
    if left >= right:
        return 0
    delta_pos = np.diff(([max(interval[0], step_pos[left])]
                            + step_pos[left + 1:right]
                            + [min(interval[1], step_pos[right])]
                            ))
    return np.dot(delta_pos, copy_numbers[left:right])

def make_cnv_plot(sample, method=""):
    """Create a genome-wide copy number plot and save data on dist."""
    logger.info(f"CNVP start")
    genome = data.ReferenceGenome()
    logger.info(f"Read positions:\n{sample.reads[:100]}")

    if method in ["", "local_coverage"]:
        bin_midpoints, cnv = mean_local_coverage(sample.reads, genome)
    elif method == "bins":
        bin_midpoints, cnv = get_cnv_from_read_start(sample.reads, genome)
    else:
        raise AttributeError(f"method has no instance {method}.")
    logger.info(f"Bin midpoints:\n{bin_midpoints}")
    logger.info(f"CNV:\n{cnv}")

    avg_read_per_bin = len(sample.reads) // len(bin_midpoints)

    cnv_plot = cnv_plot_from_data(
        data_x=bin_midpoints,
        data_y=cnv,
        avg_y = avg_read_per_bin,
        sample=sample,
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
            config.NANODIP_REPORTS,
            f"{sample_name}_{reference_name}",
        )

    def make_umap_plot(self):
        """Invoke umap plot algorithm and save to disk."""
        self.sample = data.SampleData(self.sample_name)
        self.reference = data.ReferenceData(self.reference_name) # time: 3.6s
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
        )[:config.UMAP_PLOT_TOP_MATCHES + 1]
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
            config.NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         config.ENDINGS["ranking"],
                        ),
        )
        convert_html_to_pdf(html_report, file_path)

        file_path = os.path.join(
            config.NANODIP_REPORTS,
            "%s%s" % (self.sample_name, config.ENDINGS["cpg_cnt"]),
        )
        with open(file_path, "w") as f:
            f.write("%s" % len(self.sample.cpg_overlap))

    def save_to_disk(self):
        # Save Methylation Matrix.
        file_path = os.path.join(config.NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         config.ENDINGS["methyl"])
        )
        np.save(file_path, self.methyl_overlap)

        # Save UMAP Matrix.
        file_path = os.path.join(config.NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         config.ENDINGS["umap_csv"])
        )
        self.umap_df.to_csv(file_path, index=False)

        # Write UMAP plot to disk.
        file_path = os.path.join(
            config.NANODIP_REPORTS,
            f"%s_%s%s" % (self.sample_name, self.reference_name,
                          config.ENDINGS["umap_all"]),
        )
        self.plot.write_html(file_path, config=dict({"scrollZoom": True}))
        self.plot.write_json(file_path[:-4] + "json")
        self.plot.write_image(file_path[:-4] + "png") # Time consumption 1.8s

        # Write UMAP close-up plot to disk.
        file_path = os.path.join(
            config.NANODIP_REPORTS,
            f"%s_%s%s" % (self.sample_name, self.reference_name,
                          config.ENDINGS["umap_top"]),
        )
        self.cu_plot.write_html(file_path, config=dict({"scrollZoom": True}))
        self.cu_plot.write_json(file_path[:-4] + "json")
        self.cu_plot.write_image(file_path[:-4] + "png") # Time consumption 0.9s

        # Save close up ranking report.
        self.save_ranking_report() # Time consumption 0.4s


    def files_on_disk(self):
        methyl_overlap_path = os.path.join(config.NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         config.ENDINGS["methyl"])
        )
        plot_path = self.path + config.ENDINGS["umap_all_json"]
        cu_plot_path = self.path + config.ENDINGS["umap_top_json"]

        return (os.path.exists(plot_path) and
                os.path.exists(cu_plot_path) and
                os.path.exists(methyl_overlap_path))

    def read_from_disk(self):
        methyl_overlap_path = os.path.join(config.NANODIP_REPORTS,
            "%s_%s%s" % (self.sample_name,
                         self.reference_name,
                         config.ENDINGS["methyl"])
        )
        plot_path = self.path + config.ENDINGS["umap_all_json"]
        cu_plot_path = self.path + config.ENDINGS["umap_top_json"]

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
    genome = data.ReferenceGenome()

    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.path = os.path.join(config.NANODIP_REPORTS, f"{self.sample_name}")

    def read_from_disk(self):
        plot_path = self.path + config.ENDINGS["cnv_json"]
        bins_plot_path = self.path + config.ENDINGS["cnv_bins_json"]
        genes_path = self.path + config.ENDINGS["genes"]
        with open(plot_path, "r") as f:
            self.plot_json = f.read()
        self.plot = from_json(self.plot_json)
        with open(bins_plot_path, "r") as f:
            self.bins_plot_json = f.read()
        self.genes = pd.read_csv(genes_path)
        return True

    def files_on_disk(self):
        plot_path = self.path + config.ENDINGS["cnv_json"]
        bins_plot_path = self.path + config.ENDINGS["cnv_bins_json"]
        genes_path = self.path + config.ENDINGS["genes"]
        return (
            os.path.exists(plot_path) and 
            os.path.exists(genes_path) and
            os.path.exists(bins_plot_path)
        )

    def make_cnv_plot(self):
        self.sample = data.SampleData(self.sample_name)
        self.sample.set_reads() # time consumption 2.5s
        self.reads = self.sample.reads
        self.cnv_step_pos, self.cnv_step_val = get_copy_numbers(
            CNVData.genome.length, self.reads
        )
        self.genes = self.get_gene_coverage() 
        self.relevant_genes = self.genes.loc[self.genes.relevant]
        self.plot = make_cnv_plot(self.sample) # time consumption 2.8s
        self.plot_json = self.plot.to_json()
        self.bins_plot = make_cnv_plot(self.sample, method="bins") # time consumption 2.4s
        self.bins_plot_json = self.bins_plot.to_json()
        self.save_to_disk()

    def save_to_disk(self):
        self.plot.write_html(self.path + config.ENDINGS["cnv_html"],
                             config=dict({"scrollZoom": True}))
        write_json(self.plot, self.path + config.ENDINGS["cnv_json"])
        write_json(self.bins_plot, self.path + config.ENDINGS["cnv_bins_json"])
        if not os.path.exists(self.path + config.ENDINGS["cnv_png"]):
            # time consuming operation (1.96s)
            self.plot.write_image(
                self.path + config.ENDINGS["cnv_png"], width=1280, height=720
            )
        with open(self.path + config.ENDINGS["aligned_reads"], "w") as f:
            f.write("%s" % len(self.reads))
        with open(self.path + config.ENDINGS["reads_csv"], "w") as f:
            write = csv.writer(f)
            write.writerows(self.sample.reads)
        cols = [
            "name", "seqname", "start", "end", "len", "midpoint" ,
            "cn_obs", "cn_expt", "cn_norm",
        ]
        self.genes[cols].to_csv(self.path + config.ENDINGS["genes"], index=False)
        self.relevant_genes[cols].to_csv(self.path + config.ENDINGS["relevant_genes"], index=False)

    def get_gene_coverage(self):
        genes_path = self.path + config.ENDINGS["genes"]
        #if os.path.exists(genes_path):
        #    return pd.read_csv(genes_path)
        genes = CNVData.genome.genes
        genes["interval"] = list(zip(genes.start, genes.end))
        genes["cn_obs"] = genes.interval.apply(
            lambda z: copy_number_sum(self.cnv_step_pos, self.cnv_step_val, z)
        )
        read_len_sum = sum([i[1] - i[0] for i in self.reads])
        genes["cn_expt"] = genes.interval.apply(
            lambda z: (z[1] - z[0])/CNVData.genome.length*read_len_sum
        )
        genes["cn_norm"] = genes.apply(
            lambda z: (z["cn_obs"] - z["cn_expt"]) / z["cn_expt"], axis=1)
        genes = genes.sort_values(by="cn_norm", ascending=False)
        return genes

    def get_gene_positions(self, genes):
        gene_pos = self.genes.loc[CNVData.genome.genes.name.isin(genes)]
        return gene_pos[["name", "midpoint", "cn_norm"]]

    def plot_cnv_and_genes(self, genes):
        add_scatter = self.get_gene_positions(genes)
        plot = go.Figure(self.plot)
        plot.add_trace(
            go.Scatter(
                name="",
                x=add_scatter.midpoint,
                y=add_scatter.cn_norm,
                #hovertext=genes,
                customdata=[genes],
                hovertemplate="<b>Gene: </b> %{customdata[0]} <br>",
                mode="markers",
                marker_color="rgba(100,255,0,1)",
                showlegend=False,
            ))
        return plot.to_json()


# create and cofigure logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(levelname)s %(asctime)s %(lineno)d - %(message)s")
file_handler = logging.FileHandler(
    os.path.join(config.NANODIP_REPORTS, "nanodip.log"),
    "w",
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

## output to console
# stream_handler = logging.StreamHandler()
# logger.addHandler(stream_handler)
