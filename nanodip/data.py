"""
## Data

Data containers for sample, reference-data and reference-genome/gene
data.
"""


# start_external_modules
import base64
import hashlib
import logging
import os
import re

from tqdm import tqdm
import numpy as np
import pandas as pd
import pysam
# end_external_modules

# start_internal_modules
from nanodip.config import (
    ANNOTATIONS,
    ANNOTATION_ACRONYMS_BASEL,
    ANNOTATION_ACRONYMS_TCGA,
    BETA_VALUES,
    CHROMOSOMES,
    EMPTY_SAMPLE,
    ENDING,
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
from nanodip.utils import (
    files_by_ending,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)

def binary_reference_data_exists():
    """Check if the binary form of the reference data has already been
    created.
    """
    return (
        os.path.exists(REFERENCE_METHYLATION_DATA) and
        os.path.exists(REFERENCE_METHYLATION) and
        os.path.exists(REFERENCE_CPG_SITES) and
        os.path.exists(REFERENCE_SPECIMENS) and
        os.path.exists(REFERENCE_METHYLATION_SHAPE)
    )

def make_binary_reference_data(
    input_dir=BETA_VALUES,
    output_dir=REFERENCE_METHYLATION_DATA,
    cutoff=METHYLATION_CUTOFF,
):
    """Create binary methylation files from raw reference data.

    Args:
        input_dir: Directory of raw reference data (beta values)
            as float array-files.
        output_dir: Output directory where binary methylation file
            and metadata will be written.
        cutoff: Empirical cutoff value for methylated
            (round to 1) and unmethylated (round to 0) CpGs.
    """
    print("The binary reference data is generated. Takes 5-10 minutes.")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    specimens = [f for f in os.listdir(input_dir) if f.endswith(".bin")]

    # Get shape parameters of output_data
    specimen_path0 = os.path.join(input_dir, specimens[0])
    with open(specimen_path0, "r") as f:
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
    specimen_names = [s[:-(len(ENDING["betas_bin"]) + 1)] for s in specimens]
    with open(specimens_file, "w") as f:
        f.write("\n".join(specimen_names))

    # write reference cpg sites
    index_file = os.path.join(output_dir, "cpg_sites.csv")
    with open(os.path.join(input_dir, "index.csv")) as f:
        index = f.read()
    with open(index_file, "w") as f:
        f.write(index)

def make_binary_reference_data_if_needed():
    """Creates binary reference data if not found on disk."""
    if not binary_reference_data_exists():
        make_binary_reference_data()

def _get_annotation(name, mclasses=None):
    """Reads annotation as csv file from disk, and returns it as
    pd.DataFrame. If csv is missing or file not up to date, annotation
    is read from original excel file (slow) and csv file is written to
    disk.
    """
    path_csv = os.path.join(ANNOTATIONS, name + ".csv")
    path_xlsx = os.path.join(ANNOTATIONS,name + ".xlsx")
    csv_exists_and_up_to_date = (
        os.path.exists(path_csv) and
        os.path.getmtime(path_csv) > os.path.getmtime(path_xlsx)
    )
    if csv_exists_and_up_to_date:
        annotation = pd.read_csv(path_csv)
    else:
        annotation = pd.read_excel(
            path_xlsx,
            header=None,
            names=["id", "methylation_class", "custom_text"],
            engine="openpyxl",
        )
        annotation.to_csv(path_csv, index=False)
    if mclasses is not None:
        return annotation[annotation.methylation_class.isin(mclasses)]
    return annotation

def hash_from_string(string):
    sha256_hash = hashlib.sha256(string.encode()).digest()
    filename_hash = base64.urlsafe_b64encode(sha256_hash).decode().rstrip("=")
    return filename_hash


class Reference:
    """Container of reference data and metadata."""

    # Save cpgs with index as dictionary to allow fast index lookup.
    with open(REFERENCE_CPG_SITES, "r") as f:
        cpg_site_to_index = {
            cpg:i for i, cpg in enumerate(f.read().splitlines())
        }
    # All possible CpG sites
    cpg_sites = cpg_site_to_index.keys()
    # Id's of all reference specimens
    with open(REFERENCE_SPECIMENS) as f:
        all_specimens = f.read().splitlines()
    # Save as dictionary to allow fast index lookup.
    specimen_to_index = {
        s:i for i, s in enumerate(all_specimens)
    }

    def __init__(self, name, mclasses=None):
        make_binary_reference_data_if_needed()
        self.name = name
        if mclasses is not None:
            self.name += "-" + hash_from_string("".join(mclasses))
        self.annotation = _get_annotation(name, mclasses)
        # Only consider specimens with annotation entry and binary file.
        annotated_specimens = set(self.annotation["id"]) & set(
            Reference.all_specimens
        )
        self.specimens_index = [
            Reference.specimen_to_index[a] for a in annotated_specimens
        ]
        self.specimens_index.sort()
        # Save as dictionary to allow fast methylation class lookup.
        specimen_to_mc = dict(
            zip(self.annotation.id, self.annotation.methylation_class)
        )
        # Annotated specimens sorted by increasing index
        self.specimens = [
            Reference.all_specimens[i] for i in self.specimens_index
        ]
        self.methylation_class = [
            specimen_to_mc[s] for s in self.specimens
        ]
        self.description = Reference.get_description(
            self.methylation_class
        )

    def get_description(methylation_classes):
        """Returns a description of the methylation class using
        a heuristic approach.
        """
        abbr_df = pd.read_csv(ANNOTATION_ACRONYMS_BASEL)
        abbr = dict(
            zip(abbr_df.MethylClassStr, abbr_df.MethylClassShortDescr)
        )
        non_trivial_abbr = abbr.copy()
        non_trivial_abbr.pop("-")
        tcga_df = pd.read_csv(ANNOTATION_ACRONYMS_TCGA, delimiter="\t")
        tcga = {r[1]:r[2] for r in tcga_df.itertuples()}
        def description(mc):
            """Returns description of methylation class {mc}."""
            mc = mc.upper()
            # Exact match
            if mc in abbr:
                return abbr[mc]
            # Else choose longest substring from Basel-Annotations/TCGA
            basel_substring = [a for a in non_trivial_abbr if a in mc]
            basel_substring.sort(key=len)
            tcga_substring = [a for a in tcga if a in mc]
            tcga_substring.sort(key=len)
            # Prefer Basel Annotation
            if (
                basel_substring and (
                    not tcga_substring or
                    len(basel_substring[-1]) >= len(tcga_substring[-1])
                )
            ):
                return abbr[basel_substring[-1]]
            # Else use TCGA Annotation
            if tcga_substring:
                return tcga[tcga_substring[-1]]
            # No proper annotation for "PITUI"
            if mc == "PITUI":
                return "Pituicytoma"
            return ""
        mc_description = [
            description(mc).capitalize() for mc in methylation_classes
        ]
        return mc_description

    def __str__(self):
        """Prints overview of object for debugging purposes."""
        lines = [
            f"Reference object:",
            f"name: '{self.name}'",
            f"annotation:\n{self.annotation}",
            f"specimens :\n{pd.DataFrame(self.specimens)}",
            f"specimens_index\n{pd.DataFrame(self.specimens_index)}",
            f"methylation_class: {pd.DataFrame(self.methylation_class)}",
            f"description: {pd.DataFrame(self.description)}",
            f"\n\nShared values of all Reference objects:",
            f"cpg_site_to_index:\n{pd.DataFrame(Reference.cpg_site_to_index.items())}",
            f"cpg_sites:\n{pd.DataFrame(Reference.cpg_sites)}",
            f"all_specimens:\n{pd.DataFrame(Reference.all_specimens)}",
            f"specimen_to_index:\n{pd.DataFrame(Reference.specimen_to_index.items())}",
        ]
        return "\n".join(lines)

class Genome:
    """Data container for reference genome data."""
    def __init__(self):
        self.chrom = pd.read_csv(CHROMOSOMES, delimiter="\t", index_col=False)
        self.chrom["offset"] = [0] + np.cumsum(self.chrom["len"]).tolist()[:-1]
        self.chrom["center"] = self.chrom["offset"] + self.chrom["len"]//2
        self.chrom["centromere_offset"] = (
            self.chrom["offset"]
            + (self.chrom["centromere_start"] + self.chrom["centromere_end"])
            // 2
        )
        self.length = (
            self.chrom["offset"].iloc[-1] + self.chrom["len"].iloc[-1]
        )
        if not os.path.exists(GENES):
            self.write_genes_csv()
        self.genes = pd.read_csv(GENES, delimiter="\t")

    def __iter__(self):
        """Enables looping over chromosomes."""
        return self.chrom.itertuples()

    def __len__(self):
        return self.length

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
        # Make data compatible with pythonic notation
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

    def __str__(self):
        """Prints overview of object for debugging purposes."""
        lines = [
            "Genome object:",
            f"length: {self.length}",
            f"chrom:\n{self.chrom}",
            f"genes:\n{self.genes}",
        ]
        return "\n".join(lines)

def cpg_methyl_from_reads(sample_name):
    """Returns all Illumina methylation CpG-sites with methylation
    status extracted so far.

    Args:
        sample_name: sample name to be analysed

    Returns:
        Pandas Data Frame containing the reads Illumina cpg_sites and
        methylation status.
    """
    methylation_info = pd.DataFrame(columns=["cpg_site", "methylation"])
    cpg_files = files_by_ending(
        NANODIP_OUTPUT, sample_name, ending=ENDING["methoverl_tsv"]
    )
    for f in cpg_files:
        # Some fast5 files do not contain any CpGs.
        try:
            cpgs = pd.read_csv(f, delimiter="\t", header=None,
                                names=["cpg_site", "methylation"])
            methylation_info = methylation_info.append(cpgs)
        except FileNotFoundError:
            logger.exception("Empty file encountered, skipping")
    return methylation_info.reset_index(drop=True)

class Sample:
    """Container of sample data."""
    def __init__(self, _name, cpgs=None):
        # Convert "" to EMPTY_SAMPLE
        name = EMPTY_SAMPLE if _name is "" else _name
        # Either Sample is initialized by name/id or by CpG's
        if (name is EMPTY_SAMPLE and cpgs is None) or (
            name is not EMPTY_SAMPLE and cpgs is not None
        ):
            raise ValueError("Either 'name' or 'cpgs' must be given")
        self.name = name
        self.methyl_df = None
        self.cpg_overlap = None
        self.cpg_overlap_index = None
        self.reads = None
        self.set_cpgs(cpgs)

    @classmethod
    def from_cpgs(cls, cpgs):
        """Constructor for manually define sample by CpG site set."""
        return cls(EMPTY_SAMPLE, cpgs)

    def cpgs_only(self):
        """Returns true iff only CpG set is given without methylation info."""
        return self.name is EMPTY_SAMPLE

    def set_reads(self):
        """Calculate all read start and end positions and save data
        as list to self.reads.
        """
        genome = Genome()
        bam_files = files_by_ending(NANODIP_OUTPUT, self.name, ending=".bam")
        read_positions = []
        for f in bam_files:
            samfile = pysam.AlignmentFile(f, "rb")
            # If there is no bam.bai index file, pysam will fail.
            if not samfile.has_index():
                logger.warning("No index file for %s. Skip.", f)
                continue
            for chrom in genome:
                for read in samfile.fetch(chrom.name):
                    read_positions.append([
                        # Coordinates in pysam are always 0-based (following
                        # the python convention). SAM text files use 1-based
                        # coordinates.
                        read.reference_start + chrom.offset,
                        # reference_end equals first position after alignment
                        # (following the python convention).
                        read.reference_end + chrom.offset,
                    ])
                    assert (read.reference_length != 0), "Empty read"
        self.reads = read_positions

    def set_cpgs(self, cpgs=None):
        """Sets CpG sites and additional data. This can be set manually
        (used for dummy samples containing CpG's only) or will be set
        automatically by searching the disk for methylation data
        coresponding to {sample.name}.
        """
        if cpgs is None:
            self.methyl_df = cpg_methyl_from_reads(self.name)
        else:
            self.methyl_df = pd.DataFrame([(x, None) for x in cpgs])
        self.methyl_df.columns = ["cpg_site", "methylation"]
        self.set_cpg_overlap()

    def set_cpg_overlap(self):
        """Sets CpG overlap and cpg-site-index between sample
        and reference CpGs.

        This is necessary since some probes have been skipped from the
        reference set, e.g. sex chromosomes.
        """
        self.cpg_overlap = set(self.methyl_df["cpg_site"]).intersection(
            Reference.cpg_sites)
        self.cpg_overlap_index = [
            Reference.cpg_site_to_index[f] for f in self.cpg_overlap
        ]
        self.cpg_overlap_index.sort()

    def __str__(self):
        """Prints overview of object for debugging purposes."""
        lines = [
            f"Sample object:",
            f"name: '{self.name}':",
            f"methyl_df: {self.methyl_df}",
            f"cpg_overlap: {pd.DataFrame(self.cpg_overlap)}",
            f"cpg_overlap_index: {pd.DataFrame(self.cpg_overlap_index)}",
            f"reads: {pd.DataFrame(self.reads)}",
        ]
        return "\n".join(lines)

def reference_methylation_from_index(reference_index, cpg_index):
    """Extract and return (reference-specimen x CpG-site) methylation
    submatrix from reference data.

    Args:
        reference_index: Indices of references to extract from reference
            data.
        cpg_index: Indices of Illumina CpG's to extract data.

    Returns:
        Numpy array matrix containing submatrix of reference methylation
        with rows=reference_index and columns=cpg_index.
    """
    make_binary_reference_data_if_needed()
    shape = [len(reference_index), len(cpg_index)]
    delta_offset = np.diff(reference_index, prepend=-1) - 1
    reference_submatrix = np.empty(shape, dtype=bool)
    with open(REFERENCE_METHYLATION_SHAPE, "r") as f:
        number_of_cpgs = [int(s) for s in f.read().splitlines()][1]
    with open(REFERENCE_METHYLATION, "rb") as f:
        for i, d in enumerate(delta_offset):
            reference_submatrix[i] = np.fromfile(
                f, dtype=bool, offset=d*number_of_cpgs, count=number_of_cpgs
            )[cpg_index]
    return reference_submatrix

def get_reference_methylation(reference_specimens, cpgs):
    """Extract and return (reference-specimen x CpG-site) methylation
    matrix.

        Args:
            reference_specimens: Iterable object of reference specimen id's.
            cpgs: Iterable object of CpG site id's.

        Returns:
            Binary matrix of dimension len(reference_specimens)xlen(cpgs)
            containing corresponding Methylation status.
    """
    if not cpgs:
        raise ValueError("Set of CpG sites is empty")
    if not reference_specimens:
        raise ValueError("Set of reference specimens is empty")

    reference_index = [
        Reference.specimen_to_index[a] for a in reference_specimens
    ]
    reference_index.sort()

    cpg_index = [
        Reference.cpg_site_to_index[c] for c in cpgs
    ]
    cpg_index.sort()
    return reference_methylation_from_index(reference_index, cpg_index)

def get_sample_methylation(sample):
    """Calculate and return sample methylation from reads.

    Args:
        sample: Sample to be analysed.

    Returns:
        Numpy array containing sample methylation on CpG overlap
        sites. A site is considered methylated if the mean methylation
        over all reads is greater than METHYLATION_CUTOFF.
    """
    if sample.methyl_df.empty:
        raise ValueError(
            "CpG set of {sample.name} is empty."
        )
    sample_methylation = np.full(
        len(Reference.cpg_sites), 0, dtype=bool
    )
    sample_mean_methylation = sample.methyl_df.groupby(
        "cpg_site",
        as_index=False,
    ).mean()
    for row in sample_mean_methylation.itertuples():
        cpg = row.cpg_site
        if cpg in sample.cpg_overlap:
            i = Reference.cpg_site_to_index[cpg]
            sample_methylation[i] = row.methylation > METHYLATION_CUTOFF
    return sample_methylation[sample.cpg_overlap_index]
