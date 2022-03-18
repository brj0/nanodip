import bisect
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.io import write_json, from_json
import pysam
import random
import re
import sys
import time

from config import (
    ANALYSIS_EXCLUSION_PATTERNS,
    ANNOTATIONS,
    BARCODE_NAMES,
    BROWSER_FAVICON,
    CHERRYPY_HOST,
    CHERRYPY_PORT,
    CNV_GRID,
    CNV_URL_PREFIX,
    CNV_URL_SUFFIX,
    DATA,
    DEBUG_MODE,
    ENDINGS,
    EXCLUDED_FROM_ANALYSIS,
    F5C,
    ILUMINA_CG_MAP,
    IMAGES,
    MINIMAP2,
    NANODIP_OUTPUT,
    NANODIP_REPORTS,
    NANODIP_VERSION,
    NEEDED_NUMBER_OF_BASES,
    PLOTLY_RENDER_MODE,
    READS_PER_FILE,
    READ_CPG_RSCRIPT,
    REFERENCE_GENOME_FA,
    REFERENCE_GENOME_MMI,
    RESULT_ENDINGS,
    RSCRIPT,
    SAMTOOLS,
    THIS_HOST,
    UMAP_PLOT_TOP_MATCHES,
    VERBOSITY,
)
from data import (
    get_sample_methylation,
    SampleData,
    ReferenceData,
    ReferenceGenome,
    get_reference_methylation,
)

from plots import CNVData, UMAPData
import config, data, plots, nanodip
from nanodip import methylation_caller

sample_name = "B2021_48459_20211112_BC10"
sample_name = "B2021_48700_20211112_BC11"
reference_name = "20210721_EpiDiP_anno"
reference_name = "GSE90496_IfP01"

# sample = SampleData(sample_name)
# reference = ReferenceData(reference_name)
# genome = ReferenceGenome()

# data.make_binary_reference_data()
# cnv = CNVData(sample_name)

# umapp = UMAPData(sample_name, reference_name)
# umapp.make_umap_plot()


import grpc

from nanodip import (
    mk_manager,
    minion_positions,
    minion_positions,
    run_information,
    position_from_device_id,
    getMinKnowApiStatus,
    device_status,
)

device_id = "MN26636"


# determine if anything is running and the kind of run, via set temperature


