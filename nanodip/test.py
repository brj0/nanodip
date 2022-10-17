import bisect
# import cupy
from minknow_api.tools import protocols
import numpy as np
import logging
import os
import argparse
import grpc
import pandas as pd
from tqdm import tqdm
import math
import plotly.express as px
import plotly.graph_objects as go
from plotly.io import write_json, from_json
import pysam
import random
import re
import sys
import time
from scipy.stats import binomtest
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC

sys.path.insert(0, "/applications/nanodip")

from nanodip.config import (
    ANALYSIS_EXCLUSION_PATTERNS,
    ANNOTATIONS,
    ANNOTATION_ACRONYMS_BASEL,
    ANNOTATION_ACRONYMS_TCGA,
    BARCODE_NAMES,
    BETA_VALUES,
    CHERRYPY_HOST,
    CHERRYPY_PORT,
    CHROMOSOMES,
    CNV_GRID,
    CNV_LINK,
    DATA,
    DEBUG_MODE,
    ENDING,
    UMAP_LINK,
    EPIDIP_UMAP_COORDINATE_FILES,
    EXCLUDED_FROM_ANALYSIS,
    F5C,
    GENES,
    GENES_RAW,
    ILLUMINA_CG_MAP,
    MINIMAP2,
    NANODIP_OUTPUT,
    NANODIP_REPORTS,
    NEEDED_NUMBER_OF_BASES,
    PLOTLY_RENDER_MODE,
    READS_PER_FILE,
    REFERENCE_GENOME_FA,
    REFERENCE_GENOME_MMI,
    REFERENCE_METHYLATION_SHAPE,
    RELEVANT_GENES,
    RESULT_ENDING,
    SAMTOOLS,
    THIS_HOST,
    UMAP_PLOT_TOP_MATCHES,
)
from nanodip.utils import (
    date_time_string_now,
    files_by_ending,
    discrete_colors,
    composite_path,
    bonferroni_corrected_ci,
)
from nanodip.data import (
    get_sample_methylation,
    Sample,
    Reference,
    Genome,
    get_reference_methylation,
    reference_methylation_from_index,
)
from nanodip.plots import (
    CNVData,
    UMAPData,
    pie_chart,
)
from nanodip.webui import (
    Device,
    Devices,
    minion_positions,
    run_information,
    download_epidip_data,
    device_status,
    active_run,
    number_of_called_bases,
    run_sample_id,
    start_run,
)
from nanodip.api import (
    connection_from_device_id,
    parse_args,
    is_position_selected,
    predominant_barcode,
    methylation_caller,
)
from nanodip.classifiers import (
    fit_and_evaluate_classifiers,
    training_test_data,
    evaluate_clf,
)

import nanodip.config as config
import nanodip.data as data
import nanodip.plots as plots
import nanodip.main as main
import nanodip.utils as utils
import nanodip.api as api
import nanodip.webui as webui
import nanodip.classifiers as classifiers

print("import done")

# define logger
logger = logging.getLogger(__name__)


sample_name = "B1992_24268_20211126_BC12"
reference_name = "GSE90496_IfP01"

sample = Sample(sample_name)
reference = Reference(reference_name)

analysis_dir = "/data/nanodip_output/test_06_16/FAP37158_pass_barcode06_5fa5ea82_0/"

