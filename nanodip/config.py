"""
## Configuration
Below are system-specific parameters that may or may not require adaptation.
Many variable names are self-explanatory. The key difference between
Nanopore setups are between devices provided by ONT (MinIT incl. running the
MinIT distribution on a NVIDIA Jetson developer kit such as the AGX Xavier,
GridION) and the typical Ubuntu-based MinKNOW version on x86_64 computers. The
raw data are written into a `/data` directory on ONT-based devices while they
are found in `/var/lib/minknow/data` on x86_64 installations. Make sure to
adapt your `DATA` accordingly. There are furthermore permission
issues and special folders / files in the MinKNOW data directory. These files
/ folders should be excluded from analysis through `EXCLUDED_FROM_ANALYSIS` so
that only real run folders will be parsed. Finally, the `NANODIP_OUTPUT` is the
place in which the background methylation and alignment process will place its
results by replicating the directory hierarchy of the MinKNOW data location.
It will not duplicate the data, and these data will be much smaller than raw
run data. They can be placed anywhere in the file tree, but also inside the
MinKNOW data path within a sub-folder. If the latter is the case, make sure to
apply appropriate read/write permissions. Final reports and figures generated
by NanoDiP are written into `NANODIP_REPORTS`.
"""

# start_external_modules
import os
from enum import Enum
# end_external_modules

# start_internal_modules
# end_internal_modules

NANODIP_VERSION = 31
__version__ = "31"

# Data directories for MinKNOW and NanoDiP output.
DATA = "/data"
NANODIP_OUTPUT = os.path.join(DATA, "nanodip_output")
NANODIP_REPORTS = os.path.join(DATA, "nanodip_reports")

# Reference data
REFERENCE_DATA = "/applications/reference_data"
BETA_VALUES = os.path.join(REFERENCE_DATA, "betaEPIC450Kmix_bin")
ANNOTATIONS = os.path.join(REFERENCE_DATA, "reference_annotations")
ANNOTATIONS_ABBREVIATIONS_BASEL = os.path.join(ANNOTATIONS, "mc_anno_ifp_basel.csv")
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
ANNOTATIONS_ABBREVIATIONS_TCGA = os.path.join(ANNOTATIONS, "tcga_study_abbreviations.tsv")
ILUMINA_CG_MAP = os.path.join(REFERENCE_DATA, "minimap_data/hg19_HumanMethylation450_15017482_v1-2_cgmap.tsv")
REFERENCE_METHYLATION_DATA = os.path.join(REFERENCE_DATA, "EPIC450K")
REFERENCE_METHYLATION = os.path.join(REFERENCE_METHYLATION_DATA, "methylation.bin")
REFERENCE_CPG_SITES = os.path.join(REFERENCE_METHYLATION_DATA, "cpg_sites.csv")
REFERENCE_SPECIMENS = os.path.join(REFERENCE_METHYLATION_DATA, "specimens.csv")
REFERENCE_METHYLATION_SHAPE = os.path.join(REFERENCE_METHYLATION_DATA, "shape.csv")

# Genome reference data
CHROMOSOMES = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19_chromosomes.tsv")

# Human reference genome in fa/minimap2 mmi format.
REFERENCE_GENOME_FA = "/applications/reference_data/minimap_data/hg19.fa"
REFERENCE_GENOME_MMI = "/applications/reference_data/minimap_data/hg19_20201203.mmi"

# HG19 Gene data
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
GENES_RAW = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19.refGene.gtf")
GENES = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19_genes.csv") # TODO crash if file changed
RELEVANT_GENES = os.path.join(REFERENCE_DATA, "hg19_cnv", "relevant_genes.csv")

# Beta values above cutoff will be interpreted as methylated.
METHYLATION_CUTOFF = 0.35

# Barcode strings, currently kit SQK-RBK004.
BARCODE_NAMES = [
    "barcode01","barcode02","barcode03",
    "barcode04","barcode05","barcode06",
    "barcode07","barcode08","barcode09",
    "barcode10","barcode11","barcode12",
]

# Number of reads per file. 400 works well on the Jetson AGX. Higher numbers
# increase batch size and RAM usage, lower numbers use more I/O resouces due
# to more frequent reloading of alignment reference.
READS_PER_FILE = "400"

# Number of basecalled bases until run termination occurs.
NEEDED_NUMBER_OF_BASES = 150_000_000

# URL prefix/suffix to load PDF with CNV plot for a given Sentrix ID.
CNV_URL_PREFIX = "http://s1665.rootserver.io/umapplot01/"
CNV_URL_SUFFIX = "_CNV_IFPBasel_annotations.pdf"

CNV_GRID = "/applications/reference_data/hg19_cnv/grid.json" # TODO to /tmp/nanodip

# Number of reference cases to be shown in subplot including copy
# number profile links (not advisable >200, plotly will become really
# slow)
UMAP_PLOT_TOP_MATCHES = 100

PLOTLY_RENDER_MODE = "webgl"

ANALYSIS_EXCLUSION_PATTERNS = ["_TestRun_"]
# List of files and folders in DATA to be exluded from analysis.
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
    "user_scripts",
]

# List of file name sections that identify past runs.
RESULT_ENDINGS = {
    "cnv_png": "_CNVplot.png",
    "ranking": "_NanoDiP_ranking.pdf",
    "report": "_NanoDiP_report.pdf",
    "umap_all": "_UMAP_all.html",
    "umap_top": "_UMAP_top.html",
}

ENDINGS = {
    **RESULT_ENDINGS,
    "aligned_reads": "_alignedreads.txt",
    "cnv_bins_json": "_CNV_binsplot.json",
    "cnv_html": "_CNVplot.html",
    "cnv_json": "_CNVplot.json",
    "cnv_png": "_CNVplot.png",
    "cpg_cnt":"_cpgcount.txt",
    "genes": "_genes.csv",
    "methyl": "_methyl_overlap.npy",
    "reads_csv": "_reads.csv",
    "relevant_genes": "_relevant_genes.csv",
    "umap_all_html": "_UMAP_all.html",
    "umap_all_json": "_UMAP_all.json",
    "umap_all_png": "_UMAP_all.png",
    "umap_csv": "_UMAP.csv",
    "umap_top_json": "_UMAP_top.json",
}

DEBUG_MODE = True
# 0=low log verbosity, 1=high log verbosity (with timestamps, for benchmarking and debugging)
VERBOSITY = 0 # TODO replace by logger

# Host and port on which the NanoDiP UI will be served
CHERRYPY_HOST = "localhost"
THIS_HOST = "localhost"
CHERRYPY_PORT = 8080

# The web browser favicon file for this application.
BROWSER_FAVICON = "/applications/nanodip/nanodip/static/img/favicon.ico"

# The location where image files for the web application are stored.
IMAGES ="/applications/nanodip"

# Paths to binaries for methylation calling.
F5C = "/applications/f5c/f5c"
MINIMAP2 = "/applications/nanopolish/minimap2/minimap2"
SAMTOOLS = "/applications/samtools/samtools"
# TODO del: R script that reads CpGs into simplified text file (absolute path)
RSCRIPT = "/applications/R-4.0.3/bin/Rscript"
READ_CPG_RSCRIPT="/applications/nanodip/readCpGs_mod02.R"

