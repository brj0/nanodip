import os

NANODIP_VERSION = 24

# Data directories
NANODIP_DATA = "/data"
NANODIP_OUTPUT = os.path.join(NANODIP_DATA, "nanodip_output")
NANODIP_REPORTS = os.path.join(NANODIP_DATA, "nanodip_reports")
REFERENCE_DATA = "/applications/reference_data"
BETA_VALUES = os.path.join(REFERENCE_DATA, "betaEPIC450Kmix_bin")
ANNOTATIONS = os.path.join(REFERENCE_DATA, "reference_annotations")
ANNOTATIONS_ABBREVIATIONS = "/applications/reference_data/reference_annotations/mc_anno_ifp_basel.csv"
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
TCGA_ANNOTATIONS_ABBREVIATIONS = "/applications/reference_data/reference_annotations/tcga_study_abbreviations.tsv"

# Reference data
ILUMINA_CG_MAP = os.path.join(REFERENCE_DATA, "minimap_data/hg19_HumanMethylation450_15017482_v1-2_cgmap.tsv")
REFERENCE_METHYLATION_DATA = os.path.join(REFERENCE_DATA, "EPIC450K")
REFERENCE_METHYLATION = os.path.join(REFERENCE_METHYLATION_DATA, "methylation.bin")
REFERENCE_CPG_SITES = os.path.join(REFERENCE_METHYLATION_DATA, "cpg_sites.csv")
REFERENCE_SPECIMENS = os.path.join(REFERENCE_METHYLATION_DATA, "specimens.csv")
REFERENCE_METHYLATION_SHAPE = os.path.join(REFERENCE_METHYLATION_DATA, "shape.csv")

# Genome reference data
CHROMOSOMES = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19_chromosomes.tsv")

# HG19 Gene data
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
GENES_RAW = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19.refGene.gtf")
GENES = os.path.join(REFERENCE_DATA, "hg19_cnv", "hg19_genes.csv")
RELEVANT_GENES = os.path.join(REFERENCE_DATA, "hg19_cnv", "relevant_genes.csv")

# Beta values above cutoff will be interpreted as methylated.
METHYLATION_CUTOFF = 0.35

# Number of basecalled bases until run termination occurs.
NEEDED_NUMBER_OF_BASES = 150_000_000

# URL prefix/suffix to load PDF with CNV plot for a given Sentrix ID.
CNV_URL_PREFIX = "http://s1665.rootserver.io/umapplot01/"
CNV_URL_SUFFIX = "_CNV_IFPBasel_annotations.pdf"

# Number of reference cases to be shown in subplot including copy
# number profile links (not advisable >200, plotly will become really
# slow)
UMAP_PLOT_TOP_MATCHES = 100

PLOTLY_RENDER_MODE = "webgl"

ANALYSIS_EXCLUSION_PATTERNS = ["_TestRun_"]


# List of file name sections that identify past runs.
RESULT_ENDINGS = {
    "umap_top": "_UMAP_top.html",
    "umap_all": "_UMAP_all.html",
    "report": "_NanoDiP_report.pdf",
    "cnv_png": "_CNVplot.png",
    "ranking": "_NanoDiP_ranking.pdf"
}

ENDINGS = {
    **RESULT_ENDINGS,
    "cpg_cnt":"_cpgcount.txt",
    "methyl": "_methyl_overlap.npy",
    "umap_csv": "_UMAP.csv",
    "cnv_json": "_CNVplot.json",
    "cnv_bins_json": "_CNV_binsplot.json",
    "cnv_html": "_CNVplot.html",
    "umap_all_json": "_UMAP_all.json",
    "umap_top_json": "_UMAP_top.json",
    "genes": "_genes.csv",
    "relevant_genes": "_relevant_genes.csv",
    "aligned_reads": "_alignedreads.txt", 
    "reads_csv": "_reads.csv",
}

DEBUG_MODE = True

# Host and port on which the NanoDiP UI will be served
HOST = "localhost"
PORT = 8080

# The web browser favicon file for this application.
BROWSER_FAVICON = "/applications/nanodip/favicon.ico"

# The location where image files for the web application are stored.
IMAGES ="/applications/nanodip"

# Number of reads per file. 400 works well on the Jetson AGX. Higher numbers
# increase batch size and RAM usage, lower numbers use more I/O resouces due
# to more frequent reloading of alignment reference.
READS_PER_FILE = "400"
