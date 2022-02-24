from tqdm import tqdm
import csv
import numpy as np
import openpyxl
import os
import pandas as pd
import pysam
import re
import time

import config

def binary_reference_data_exists():
    """Check if the binary form of the reference data was already created."""
    return (
        os.path.exists(config.REFERENCE_METHYLATION_DATA) and
        os.path.exists(config.REFERENCE_METHYLATION) and
        os.path.exists(config.REFERENCE_CPG_SITES) and
        os.path.exists(config.REFERENCE_SPECIMENS) and
        os.path.exists(config.REFERENCE_METHYLATION_SHAPE)
    )

def make_binary_reference_data(input_dir=config.BETA_VALUES,
                               output_dir=config.REFERENCE_METHYLATION_DATA,
                               cutoff=config.METHYLATION_CUTOFF):
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
        annotation_path = os.path.join(config.ANNOTATIONS, name + ".xlsx")
        self.annotation = pd.DataFrame(
            openpyxl.load_workbook(annotation_path).active.values,
            columns=["id", "methylation_class", "custom_text"],
        )
        with open(config.REFERENCE_CPG_SITES, "r") as f:
        # save as Dictionary to allow fast index lookup
            self.cpg_sites = {cpg:i for i, cpg in enumerate(
                f.read().splitlines()
            )}

        with open(config.REFERENCE_SPECIMENS) as f:
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
        methyl_dict = {row["id"]:row["methylation_class"]
            for _, row in self.annotation.iterrows()}
        self.specimen_ids = [self.specimens[i]
            for i in self.annotated_specimens_index]
        self.methylation_class = [methyl_dict[s] for s in self.specimen_ids]

        abbr_df = pd.read_csv(config.ANNOTATIONS_ABBREVIATIONS)
        abbr = {row["MethylClassStr"]:row["MethylClassShortDescr"]
                     for _, row in abbr_df.iterrows()} # TODO
        non_trivial_abbr = abbr.copy()
        non_trivial_abbr.pop("-")
        tcga_df = pd.read_csv(config.TCGA_ANNOTATIONS_ABBREVIATIONS, delimiter="\t")
        tcga = {r[0]:r[1] for _, r in tcga_df.iterrows()}
        def description(mc):
            if mc in abbr:
                return abbr[mc]
            substrings = [a for a in non_trivial_abbr if a in mc]
            substrings.sort(key=lambda x: len(x))
            if substrings:
                return abbr[substrings[-1]]
            tcga_substring = [a for a in tcga if a in mc]
            tcga_substring.sort(key=lambda x: len(x))
            if tcga_substring:
                return tcga[tcga_substring[-1]]
            else:
                return ""
        self.mc_description = [description(mc).capitalize() for mc in self.methylation_class]
        #TODO BUG ('TCGA-THCATCGA_THCA_CLASS', 'Hepatocellular adenoma') # DONE

class ReferenceGenome:

    def __init__(self):
        self.chrom = pd.read_csv(config.CHROMOSOMES,
                                 delimiter="\t",
                                 index_col=False)
        self.chrom["offset"] = [0] + np.cumsum(self.chrom["len"]).tolist()[:-1]
        self.chrom["center"] = self.chrom["offset"] + self.chrom["len"]//2
        self.chrom["centromere_offset"] = (self.chrom["offset"]
            + (self.chrom["centromere_start"] + self.chrom["centromere_end"])//2)
        self.length = (self.chrom["offset"].iloc[-1]
                     + self.chrom["len"].iloc[-1])
        if not os.path.exists(config.GENES):
            self.write_genes_csv()
        self.set_genes()

    def __iter__(self):
        return self.chrom.itertuples()

    def set_genes(self):
        """Read and set genes from csv file."""
        self.genes = pd.read_csv(
            config.GENES,
            delimiter="\t",
        )

    def write_genes_csv(self):
        """Write csv gene list with one selected transcript per gene."""
        genes = pd.read_csv(
            config.GENES_RAW,
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
        genes = genes.drop_duplicates(subset=["name", "seqname"], keep="first")
        genes = genes.sort_values("name")
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
        with open(config.RELEVANT_GENES, "r") as f:
            relevant_genes = f.read().splitlines()
        genes["relevant"] = genes.name.apply(lambda x: x in relevant_genes)
        genes["len"] = genes["end"] - genes["start"]
        genes[["name", "seqname", "start", "end", "len", "midpoint", "relevant"]].to_csv(
            config.GENES, index=False, sep="\t")

class SampleData:
    """Container of sample data"""
    def __init__(self, name):
        self.name = name
        self.cpg_sites = SampleData.get_read_cpgs(name)
        self.cpg_overlap = None
        self.cpg_overlap_index = None
        self.reads = None

    def set_reads(self, file_name=None):
        """Calculate read positions"""
        if file_name is None:
            self.reads = self.get_read_positions(ReferenceGenome())
        else:
            with open(file_name, "r") as f:
                reads = csv.reader(f)
                self.reads = [[int(r[0]),int(r[1])] for r in reads]


    def get_read_cpgs(name):
        """Get all Ilumina methylation sites with methylaton status
        within a samples reads.

        Args:
            name: sample name to be analysed

        Returns:
            Pandas Data Frame containing the reads Ilumina cpg_sites and
            methylation status.
        """

        sample_path = os.path.join(config.NANODIP_OUTPUT, name)

        if not os.path.exists(sample_path):
            raise FileNotFoundError(sample_path)

        cpg_files = []
        for root, _, files in os.walk(sample_path):
            cpg_files.extend(
                [os.path.join(root, f)
                for f in files if f.endswith("methoverlap.tsv")]
            )

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

    def get_read_positions(self, genome):
        """Return list containing start and end positions of all reads."""
        bam_files = []
        sample_path = os.path.join(config.NANODIP_OUTPUT, self.name)
        for root, _, files in os.walk(sample_path):
            bam_files.extend(
                [os.path.join(root, f)
                for f in files if f.endswith(".bam")]
            )
        read_positions = []
        for f in bam_files:
            samfile = pysam.AlignmentFile(f, "rb")
            for chrom in genome:
                for read in samfile.fetch(chrom.name):
                    read_positions.append([
                        read.reference_start + chrom.offset,
                        # reference_end equals first position after alignment
                        read.reference_end + chrom.offset,
                    ])
                    assert (read.reference_length != 0), "Empty read"
        return read_positions

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

    with open(config.REFERENCE_METHYLATION_SHAPE, "r") as f:
        number_of_cpgs = [int(s) for s in f.read().splitlines()][1]

    with open(config.REFERENCE_METHYLATION, "rb") as f:
        for i, d in enumerate(delta_offset):
            reference_matrix[i] = np.fromfile(
                f, dtype=bool, offset=d*number_of_cpgs, count=number_of_cpgs
            )[cpg_index]
    return reference_matrix

def get_reference_methylation(sample, reference):
    """Extract and return methylation information matrix from overlap of sample
    CpG's with annotated reference data."""

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
                                    config.METHYLATION_CUTOFF
    sample_methylation = sample_methylation[sample.cpg_overlap_index]
    return sample_methylation
