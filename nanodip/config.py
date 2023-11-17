"""
Configuration
===============================================================================
Below are system-specific parameters that may or may not require adaptation.
Many variable names are self-explanatory. The key difference between Nanopore
setups are between devices provided by ONT (MinIT inclusively running the MinIT
distribution on a NVIDIA Jetson developer kit such as the AGX Xavier, GridION)
and the typical Ubuntu-based MinKNOW version on x86_64 computers. The raw data
are written into a `/data` directory on ONT-based devices while they are found
in `/var/lib/minknow/data` on x86_64 installations. Make sure to adapt your
`DATA` accordingly. There are furthermore permission issues and special folders
/ files in the MinKNOW data directory. These files / folders should be excluded
from analysis through `EXCLUDED_FROM_ANALYSIS` so that only real run folders
will be parsed. Finally, the `NANODIP_OUTPUT` is the place in which the
background methylation and alignment process will place its results by
replicating the directory hierarchy of the MinKNOW data location.  It will not
duplicate the data, and these data will be much smaller than raw run data. They
can be placed anywhere in the file tree, but also inside the MinKNOW data path
within a sub-folder. If the latter is the case, make sure to apply appropriate
read/write permissions. Final reports and figures generated by NanoDiP are
written into `NANODIP_REPORTS`.
"""

# start_external_modules
import os
# end_external_modules

# start_internal_modules
# end_internal_modules


"""
General
-------------------------------------------------------------------------------
"""

# NanoDiP version number.
__version__ = "0.0.0"

# Enables/disables debug mode.
DEBUG_MODE = True


"""
Data directories for MinKNOW and NanoDiP output.
-------------------------------------------------------------------------------
"""

# Where MinKNOW places its data.
DATA = "/data"

# Location to write intermediate analysis data, i.e. methylation and alignment
# files.
NANODIP_OUTPUT = os.path.join(DATA, "nanodip_output")

# Location to write reports and figures.
NANODIP_REPORTS = os.path.join(DATA, "nanodip_reports")

# Where to write temporary files
TMP = os.path.join(DATA, "tmp")


"""
Reference data
-------------------------------------------------------------------------------
"""

# Location of all reference data.
REFERENCE_DATA = "/applications/reference_data"

# Location of preprocessed beta values.
BETA_VALUES = os.path.join(REFERENCE_DATA, "betaEPIC450Kmix_bin")

# Location of annotation spreadsheets.
ANNOTATIONS = os.path.join(REFERENCE_DATA, "reference_annotations")

# Location of annotation acronyms.
ANNOTATION_ACRONYMS = os.path.join(ANNOTATIONS, "acronyms")

# Spreadsheet containing description of annotation codes (Basel internal).
ANNOTATION_ACRONYMS_BASEL = os.path.join(
    ANNOTATION_ACRONYMS, "mc_anno_ifp_basel.csv"
)

# Spreadsheet containing description of annotation codes (TCGA). From:
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
ANNOTATION_ACRONYMS_TCGA = os.path.join(
    ANNOTATION_ACRONYMS, "tcga_study_abbreviations.tsv"
)

# Illumina probe names of the 450K array.
ILLUMINA_CPG_MAP = os.path.join(
    REFERENCE_DATA,
    "minimap_data/hg19_HumanMethylation450_15017482_v1-2_cgmap.tsv",
)

# Location of binary methylation data and metadata.
REFERENCE_METHYLATION_DATA = os.path.join(REFERENCE_DATA, "EPIC450K")

# Binary methylation file.
REFERENCE_METHYLATION = os.path.join(
    REFERENCE_METHYLATION_DATA, "methylation.bin"
)

# Metadata containing list of CpG names.
REFERENCE_CPG_SITES = os.path.join(REFERENCE_METHYLATION_DATA, "cpg_sites.csv")

# Metadata containing names of specimens.
REFERENCE_SPECIMENS = os.path.join(REFERENCE_METHYLATION_DATA, "specimens.csv")

# Metadata containing methylation matrix dimensions.
REFERENCE_METHYLATION_SHAPE = os.path.join(
    REFERENCE_METHYLATION_DATA, "shape.csv"
)

# Genome reference data containing chromosome lengths and centromere position.
CHROMOSOMES = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19_chromosomes.tsv")

# Human reference genome in fa format.
REFERENCE_GENOME_FA = os.path.join(REFERENCE_DATA, "minimap_data/hg19.fa")

# Human reference genome in minimap2 mmi format.
REFERENCE_GENOME_MMI = os.path.join(
    REFERENCE_DATA, "minimap_data/hg19_20201203.mmi"
)

# HG19 Gene data downloaded from:
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
GENES_RAW = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19.refGene.gtf")

# Contains the data from GENES_RAW in simplified form.
GENES = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19_genes.csv")

# List of clinically important genes.
RELEVANT_GENES = os.path.join(REFERENCE_DATA, "hg19_cnv", "relevant_genes.csv")


"""
File handling
-------------------------------------------------------------------------------
"""

# String patterns in sample names that exclude data from downstream analysis,
# e.g., test runs
ANALYSIS_EXCLUSION_PATTERNS = ["_TestRun_"]

# List of files and folders in DATA to be excluded from analysis.
EXCLUDED_FROM_ANALYSIS = [
    ".Trash-1000",
    "core-dump-db",
    "intermediate",
    "lost+found",
    "minimap_data",
    "nanodip_output",
    "nanodip_reports",
    "nanodip_tmp",
    "non-ont",
    "pings",
    "playback_raw_runs",
    "queued_reads",
    "raw_for_playback",
    "reads",
    "tmp",
    "user_scripts",
]

# List of file name sections that identify past runs.
RESULT_ENDING = {
    "clf_txt": "classifiers.txt",
    "cnv_png": "CNVplot.png",
    "ranking_pdf": "NanoDiP_ranking.pdf",
    "report_pdf": "NanoDiP_report.pdf",
    "umaptop_html": "UMAP_top.html",
    "umapall_html": "UMAP_all.html",
}

# Used file endings.
ENDING = {
    **RESULT_ENDING,
    "alignedreads_txt": "alignedreads.txt",
    "betas_bin": "betas_filtered.bin",
    "binmdpnts_npy": "binmidpoints.npy",
    "cnv_html": "CNVplot.html",
    "cnv_json": "CNVplot.json",
    "cnv_npy": "cnv.npy",
    "cnv_pdf": "CNVplot.pdf",
    "cpg_cnt": "cpgcount.txt",
    "freq_tsv": "-freq.tsv",
    "genes_csv": "genes.csv",
    "methoverl_npy": "methyl_overlap.npy",
    "methoverl_tsv": "-methoverlap.tsv",
    "methoverlcnt_txt": "-methoverlapcount.txt",
    "pie_png": "pie.png",
    "reads_csv": "reads.csv",
    "readsort_bam": "-reads_sorted.bam",
    "readsort_bai": "-reads_sorted.bam.bai",
    "relgenes_csv": "relevant_genes.csv",
    "result_tsv": "-result.tsv",
    "stdarr_bin": "stdarr.bin",
    "stdsortarr_bin": "stdsortarr.bin",
    "umap_csv": "UMAP.csv",
    "umap_xlsx": "UMAP.xlsx",
    "umapall_json": "UMAP_all.json",
    "umapall_png": "UMAP_all.png",
    "umaptop_json": "UMAP_top.json",
    "umaptop_png": "UMAP_top.png",
}

# Name for Samples cpgs containing only
EMPTY_SAMPLE = "nosample"

# Temporary file directory containing status information about samples
SMP_STATUS_TMP = os.path.join(TMP, "samples", "status")
SMP_MCACTIVE_TMP = os.path.join(TMP, "samples", "mcactive")


"""
Experiments & basecalling
-------------------------------------------------------------------------------
"""

# Beta values above this cutoff will be interpreted as methylated.
METHYLATION_CUTOFF = 0.35

# Barcode strings, currently kit SQK-RBK004.
BARCODE_NAMES = [
    "barcode01",
    "barcode02",
    "barcode03",
    "barcode04",
    "barcode05",
    "barcode06",
    "barcode07",
    "barcode08",
    "barcode09",
    "barcode10",
    "barcode11",
    "barcode12",
]

# Number of reads per file. 400 works well on the Jetson AGX. Higher numbers
# increase batch size and RAM usage, lower numbers use more I/O resources due
# to more frequent reloading of alignment reference.
READS_PER_FILE = "400"

# Number of basecalled bases until run termination occurs.
NEEDED_NUMBER_OF_BASES = 150_000_000

# Paths to binaries for methylation calling.
F5C = "/applications/f5c/f5c"
MINIMAP2 = "/applications/nanopolish/minimap2/minimap2"
SAMTOOLS = "/applications/samtools/samtools"


"""
UMAP/CNV plots, Epidip
-------------------------------------------------------------------------------
"""

# URL to load PDF with CNV plot for a given Sentrix ID (substituted for '%s')
CNV_LINK = (
    "http://s1665.rootserver.io/umapplot01/%s_CNV_IFPBasel_annotations.pdf"
)

# URL to load precalculated UMAP coordinates.
UMAP_LINK = "http://s1665.rootserver.io/umap_links/%s"

# Path to precalculated CNV plotly grid.
CNV_GRID = "/applications/reference_data/hg19_cnv/grid.json"

# Number of reference cases to be shown in subplot including copy
# number profile links (not advisable >200, plotly will become really
# slow)
UMAP_PLOT_TOP_MATCHES = 100

# Controls the browser API used to draw marks. Default "webgl", alternative
# "svg" without proper webgl support (e.g. for Firefox, use "svg"; slower,
# but does not require GPU).
PLOTLY_RENDER_MODE = "webgl"

# List of UMAP coordinate files hosted on EpiDiP.
EPIDIP_UMAP_COORDINATE_FILES = [
    "UMAP_all_bVals_top_25000.xlsx",
    "UMAP_all_bVals_top_50000.xlsx",
    "UMAP_all_bVals_top_75000.xlsx",
    "gpumap_25000.xlsx",
    "gpumap_50000.xlsx",
    "gpumap_75000.xlsx",
]

# EpiDiP temporary file directory
EPIDIP_TMP = os.path.join(TMP, "epidip")

# Desired RAM usage. 4 works best on Jetson AGX 32GB; adapt to GPU / RAM
# layout, see
# https://forums.developer.nvidia.com/t/nvmapreserveop-0x80000000-failed-22-when-running-cufft-plans/168781/14
GPU_RAM_USAGE = 4 * 1024**3

# Size per float in GPU RAM (tested on AGX Xavier)
GPU_FLOAT_SIZE = 8


"""
CherryPy
-------------------------------------------------------------------------------
"""

# Host and port on which the NanoDiP UI will be served
CHERRYPY_HOST = "localhost"
THIS_HOST = "localhost"
CHERRYPY_PORT = 8080
