"""
## Data

Data containers for sample, reference-data and reference-genome/gene
data.
"""


# start_external_modules
from tqdm import tqdm
import numpy as np
import os
import pandas as pd
import pysam
import re
import time
# end_external_modules

# start_internal_modules
from config import (
    ANNOTATIONS,
    ANNOTATIONS_ABBREVIATIONS_BASEL,
    ANNOTATIONS_ABBREVIATIONS_TCGA,
    BETA_VALUES,
    CHROMOSOMES,
    GENES,
    GENES_RAW,
    METHYLATION_CUTOFF,
    NANODIP_OUTPUT,
    REFERENCE_CPG_SITES,
    REFERENCE_METHYLATION,
    REFERENCE_METHYLATION_DATA,
    REFERENCE_METHYLATION_SHAPE,
    REFERENCE_SPECIMENS,
    RELEVANT_GENES,
)
from utils import (
    date_time_string_now
)
# end_internal_modules

def binary_reference_data_exists():
    """Check if the binary form of the reference data was already created."""
    return (
        os.path.exists(REFERENCE_METHYLATION_DATA) and
        os.path.exists(REFERENCE_METHYLATION) and
        os.path.exists(REFERENCE_CPG_SITES) and
        os.path.exists(REFERENCE_SPECIMENS) and
        os.path.exists(REFERENCE_METHYLATION_SHAPE)
    )

def make_binary_reference_data(input_dir=BETA_VALUES,
                               output_dir=REFERENCE_METHYLATION_DATA,
                               cutoff=METHYLATION_CUTOFF):
    """Create binary methylation files from raw reference data.

    Args:
        input_dir: Directory of reference data as float arrays-
            files.
        output_dir: Output dir containing binary array-file.
        cutoff: Empirical cutoff value for methylated
            (round to 1) and unmethylated (round to 0) CpGs.
    """
    print("The binary reference data is generated. Takes 5-10 minutes.")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    specimens = [f for f in os.listdir(input_dir)
                 if f.endswith(".bin")]

    # Get shape parameters of output_data
    path0 = os.path.join(input_dir, specimens[0])
    with open(path0, "r") as f:
        beta_values_0 = np.fromfile(f, dtype=float)
    shape = (len(specimens), len(beta_values_0))

    methylation_data = np.empty(shape, dtype=bool)

    for i, specimen in enumerate(tqdm(specimens, desc="Reading reference")):
        specimen_path = os.path.join(input_dir, specimen)

        with open(specimen_path, "rb") as f:
            beta_values = np.fromfile(f, dtype=float)
            methylation_data[i] = np.digitize(
                beta_values,
                bins=[cutoff]
            ).astype(bool)

    # write methylation data as binary
    methylation_file = os.path.join(output_dir, "methylation.bin")
    methylation_data.tofile(methylation_file)

    # write shape parameters
    shape_file = os.path.join(output_dir, "shape.csv")
    with open(shape_file, "w") as f:
        f.write("%s\n %s" % shape)

    # write reference specimens
    specimens_file = os.path.join(output_dir, "specimens.csv")
    specimen_names = [s[:-len("_betas_filtered.bin")] for s in specimens]
    with open(specimens_file, "w") as f:
        f.write("\n".join(specimen_names))

    # write reference cpg sites
    index_file = os.path.join(output_dir, "cpg_sites.csv")
    with open(os.path.join(input_dir, "index.csv")) as f:
        index = f.read()
    with open(index_file, "w") as f:
        f.write(index)

def make_binary_reference_data_if_needed():
    if not binary_reference_data_exists():
        make_binary_reference_data()

class ReferenceData:
    """Container of reference data and metadata."""
    def __init__(self, name):
        make_binary_reference_data_if_needed()
        self.name = name
        self.annotation = self.get_annotation()
        with open(REFERENCE_CPG_SITES, "r") as f:
        # save as Dictionary to allow fast index lookup
            self.cpg_sites = {cpg:i for i, cpg in enumerate(
                f.read().splitlines()
            )}

        with open(REFERENCE_SPECIMENS) as f:
            self.specimens = f.read().splitlines()

        # determine if there are entries in the annotation without corresponding
        # methylation binary file
        self.annotated_specimens = list(
            set(self.annotation["id"]) & set(self.specimens)
        )

        # Save as dictionary to allow fast hash lookup.
        index = {s:i for i, s in enumerate(self.specimens)}
        self.annotated_specimens_index = [index[a]
            for a in self.annotated_specimens]
        self.annotated_specimens_index.sort()

        # Save as dictionary to allow fast hash lookup.
        methyl_dict = {i:mc for i, mc in
            zip(self.annotation.id, self.annotation.methylation_class)
        }
        self.specimen_ids = [self.specimens[i]
            for i in self.annotated_specimens_index]
        self.methylation_class = [methyl_dict[s] for s in self.specimen_ids]
        self.description = ReferenceData.get_description(
            self.methylation_class
        )

    def get_description(methylation_classes):
        """Returns a description of the methylation class."""
        abbr_df = pd.read_csv(ANNOTATIONS_ABBREVIATIONS_BASEL)
        abbr = {
            mc:desc for mc, desc in
            zip(abbr_df.MethylClassStr, abbr_df.MethylClassShortDescr)
        }
        non_trivial_abbr = abbr.copy()
        non_trivial_abbr.pop("-")
        tcga_df = pd.read_csv(ANNOTATIONS_ABBREVIATIONS_TCGA, delimiter="\t")
        tcga = {r[0]:r[1] for _, r in tcga_df.iterrows()}
        def description(mc):
            mc = mc.upper()
            # Exact match
            if mc in abbr:
                return abbr[mc]
            # Else choose longest substring from Basel-Annotations/TCGA
            basel_substring = [a for a in non_trivial_abbr if a in mc]
            basel_substring.sort(key=lambda x: len(x))
            tcga_substring = [a for a in tcga if a in mc]
            tcga_substring.sort(key=lambda x: len(x))
            # Prefer Basel Annotation
            if (
                basel_substring and (
                    not tcga_substring or
                    len(basel_substring[-1]) >= len(tcga_substring[-1])
                )
            ):
                return abbr[basel_substring[-1]]
            if tcga_substring:
                return tcga[tcga_substring[-1]]
            # No proper annotation for "PITUI"
            if mc == "PITUI":
                return "Pituicytoma"
            else:
                return ""
        mc_description = [
            description(mc).capitalize() for mc in methylation_classes
        ]
        return mc_description

    def get_annotation(self):
        """Reads annotation as csv file from disk, and returns is as
        pd.DataFrame. If csv is missing or file not up to date, annotation
        is read from original excel file (slow) and csv file is written to
        disk.
        """
        path_csv = os.path.join(ANNOTATIONS, self.name + ".csv")
        path_xlsx = os.path.join(ANNOTATIONS, self.name + ".xlsx")
        csv_exists_and_up_to_date = (
            os.path.exists(path_csv) and
            os.path.getmtime(path_csv) > os.path.getmtime(path_xlsx)
        )
        if csv_exists_and_up_to_date:
            return pd.read_csv(path_csv)
        annotation = pd.read_excel(
            path_xlsx,
            header=None,
            names=["id", "methylation_class", "custom_text"],
            engine="openpyxl",
        )
        annotation.to_csv(path_csv, index=False)
        return annotation

class ReferenceGenome:

    def __init__(self):
        self.chrom = pd.read_csv(CHROMOSOMES,
                                 delimiter="\t",
                                 index_col=False)
        self.chrom["offset"] = [0] + np.cumsum(self.chrom["len"]).tolist()[:-1]
        self.chrom["center"] = self.chrom["offset"] + self.chrom["len"]//2
        self.chrom["centromere_offset"] = (self.chrom["offset"]
            + (self.chrom["centromere_start"] + self.chrom["centromere_end"])//2)
        self.length = (self.chrom["offset"].iloc[-1]
                     + self.chrom["len"].iloc[-1])
        if not os.path.exists(GENES):
            self.write_genes_csv()
        self.set_genes()

    def __iter__(self):
        return self.chrom.itertuples()

    def set_genes(self):
        """Read and set genes from csv file."""
        self.genes = pd.read_csv(
            GENES,
            delimiter="\t",
        )

    def write_genes_csv(self):
        """Write csv gene list with one selected transcript per gene."""
        genes = pd.read_csv(
            GENES_RAW,
            delimiter="\t",
            names=["seqname", "source", "feature", "start", "end",
                   "score", "strand", "frame", "attribute"],
            usecols=["seqname", "feature", "start", "end", "attribute"]
        )
        genes = genes.loc[
            (genes["feature"] == "transcript")
            & (genes["seqname"].isin(self.chrom.name))
        ]
        genes["name"] = genes.attribute.apply(
            lambda x: re.search('gene_name(.*)"(.*)"', x).group(2)
        )
        genes["transcript"] = genes.attribute.apply(
            lambda x: re.search(
                'transcript_id(.*)"(.*)"(.*)gene_name(.*)', x
                ).group(2)
        )
        genes = genes.drop_duplicates(subset=["name", "seqname"], keep="first")
        genes = genes.sort_values("name")
        genes["loc"] = genes.apply(
            lambda z: (
                  z["seqname"]
                + ":"
                + "{:,}".format(z["start"])
                + "-"
                + "{:,}".format(z["end"])
            ),
            axis=1,
        )
        # Make data comapitle with pythonic notation
        genes["end"] += 1
        offset = {i.name:i.offset for i in self}
        genes["start"] = genes.apply(
            lambda z: offset[z["seqname"]] + z["start"],
            axis=1,
        )
        genes["end"] = genes.apply(
            lambda z: offset[z["seqname"]] + z["end"],
            axis=1,
        )
        genes["midpoint"] = (genes["start"] + genes["end"]) // 2
        with open(RELEVANT_GENES, "r") as f:
            relevant_genes = f.read().splitlines()
        genes["relevant"] = genes.name.apply(lambda x: x in relevant_genes)
        genes["len"] = genes["end"] - genes["start"]
        genes[["name", "seqname", "start", "end",
               "len", "midpoint", "relevant", "transcript",
               "loc",
        ]].to_csv(GENES, index=False, sep="\t")

def files_by_ending(directory, sample_name, ending):
    """Returns a list containing all sample output files with a given
    ending.
    """
    sample_path = os.path.join(directory, sample_name)
    output_files = []
    for root, _, files in os.walk(sample_path):
        output_files.extend(
            [os.path.join(root, f)
            for f in files if f.endswith(ending)]
        )
    return output_files

class SampleData:
    """Container of sample data."""
    def __init__(self, name):
        self.name = name
        self.cpg_sites = SampleData.get_read_cpgs(name)
        self.cpg_overlap = None
        self.cpg_overlap_index = None
        self.reads = None

    def set_reads(self):
        """Calculate all read start and end positions."""
        genome = ReferenceGenome()
        bam_files = files_by_ending(NANODIP_OUTPUT, self.name, ending=".bam")
        read_positions = []
        for f in bam_files:
            samfile = pysam.AlignmentFile(f, "rb")
            for chrom in genome:
                for read in samfile.fetch(chrom.name):
                    read_positions.append([
                        read.reference_start + chrom.offset,
                        # reference_end equals first position after alignment
                        # consistent with python notations.
                        read.reference_end + chrom.offset,
                    ])
                    assert (read.reference_length != 0), "Empty read"
        self.reads = read_positions

    def get_read_cpgs(sample_name):
        """Get all Ilumina methylation sites with methylaton status
        within a samples reads.

        Args:
            sample_name: sample name to be analysed

        Returns:
            Pandas Data Frame containing the reads Ilumina cpg_sites and
            methylation status.
        """

        sample_path = os.path.join(NANODIP_OUTPUT, sample_name)

        if not os.path.exists(sample_path):
            raise FileNotFoundError(sample_path)

        cpg_files = files_by_ending(NANODIP_OUTPUT, sample_name,
                                     ending="methoverlap.tsv")

        methylation_info = pd.DataFrame()

        for f in cpg_files:
            # Some fast5 files do not contain any CpGs.
            try:
                cpgs = pd.read_csv(f, delimiter="\t", header=None,
                                   names=["cpg_site", "methylation"])
                methylation_info = methylation_info.append(cpgs)
            except FileNotFoundError:
                logger.exception("empty file encountered, skipping")

        return methylation_info

    def set_cpg_overlap(self, reference):
        """Calculate CpG overlap data between sample and reference.

        Some probes have been skipped from the reference set, e.g. sex
        chromosomes.
        """ #TODO is this true?
        self.cpg_overlap = set(self.cpg_sites["cpg_site"]).intersection(
            reference.cpg_sites.keys())

        self.cpg_overlap_index = [reference.cpg_sites[f]
            for f in self.cpg_overlap]
        self.cpg_overlap_index.sort()

def _get_reference_methylation(reference_index, cpg_index):
    """Extract and return methylation information matrix from reference data.

    Args:
        reference_index: Index of references to extract from reference
            data.
        cpg_index: Index of CpG's to extract from CpG data.

    Returns:
        Numpy array matrix containing submatrix of reference data
        with rows=reference_index and columns=cpg_index.
    """

    make_binary_reference_data_if_needed()
    shape = [len(reference_index), len(cpg_index)]
    delta_offset = np.diff(reference_index, prepend=-1) - 1
    reference_matrix = np.empty(shape, dtype=bool)

    with open(REFERENCE_METHYLATION_SHAPE, "r") as f:
        number_of_cpgs = [int(s) for s in f.read().splitlines()][1]

    with open(REFERENCE_METHYLATION, "rb") as f:
        for i, d in enumerate(delta_offset):
            reference_matrix[i] = np.fromfile(
                f, dtype=bool, offset=d*number_of_cpgs, count=number_of_cpgs
            )[cpg_index]
    return reference_matrix

def get_reference_methylation(sample, reference):
    """Extract and return methylation information matrix from overlap of sample
    CpG's with annotated reference data.
    """

    reference_index = reference.annotated_specimens_index
    cpg_index = sample.cpg_overlap_index
    result = _get_reference_methylation(reference_index, cpg_index)
    return result


def get_sample_methylation(sample, reference):
    """Calculate sample methylation info from reads.

    Args:
        sample: Sample data set.
        reference: Reference data set.
        cpg_ovelrap: Set containing intersection of cpg sites in cpg_sample
            and reference_cpg_site.
    Returns:
        Numpy array containing sample Methylation information.
    """

    sample_methylation = np.full(len(reference.cpg_sites), 0, dtype=bool)
    sample_mean_methylation = sample.cpg_sites.groupby(
        "cpg_site",
        as_index=False).mean()

    for _, row in sample_mean_methylation.iterrows():
        cpg = row["cpg_site"]
        if cpg in sample.cpg_overlap:
            i = reference.cpg_sites[cpg]
            sample_methylation[i] = row["methylation"] > \
                                    METHYLATION_CUTOFF
    sample_methylation = sample_methylation[sample.cpg_overlap_index]
    return sample_methylation
