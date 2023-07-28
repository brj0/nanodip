import bisect

# import cupy
from minknow_api.tools import protocols
import threading
import numpy as np
import logging
import os
import argparse
import grpc
import pandas as pd
import datetime
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

pdp = lambda x: print(x.to_string())

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
    EPIDIP_TMP,
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
from nanodip.epidip import (
    calculate_std,
    gpu_enabled,
    top_variable_cpgs,
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

# sample_name = "test20221124a"
# sample_name = "B2022_30785_20220715_BC12"
# sample = Sample(sample_name)

# reference_name = "MNG_IfP_v1"
# reference = Reference(reference_name)
# calculate_std(reference)

# reference_name = "AllIDATv2_20210804"
# reference = Reference(reference_name)
# calculate_std(reference)

reference_id = "GSE90496_IfP01"
reference = Reference(reference_id)
std_all = calculate_std(reference)


gbm_all = [
    "GBM_G34",
    "GBM_LOW",
    "GBM_MES",
    "GBM_MID",
    "GBM_MYCN",
    "GBM_NOS",
    "GBM_RTK_I",
    "GBM_RTK_II",
    "GBM_RTK_III",
]
mng_ben = ["MNG_BEN-1", "MNG_BEN-2"]
pitad_all = [
    "PITAD",
    "PITAD_FSH_LH",
    "PITAD_ACTH",
    "PITAD_STH_SPA",
    "PITAD_STH_DNS_A",
    "PITAD_TSH",
    "PITAD_STH_DNS_B",
    "PITAD_PRL",
]


brain = Reference(reference_id)
gbm = Reference(reference_id, mclasses=gbm_all)
mng = Reference(reference_id, mclasses=mng_ben)
pitad = Reference(reference_id, mclasses=pitad_all)
mng_all = Reference("MNG_IfP_v1")
mng = Reference("MNG_IfP_v1", mng_ben)


mgbm = get_reference_methylation(
    gbm.specimens, gbm.cpg_sites
)
mmng = get_reference_methylation(
    mng.specimens, mng.cpg_sites
)
mpitad = get_reference_methylation(
    pitad.specimens, pitad.cpg_sites
)

sgbm = mgbm[-1, :]
smng = mmng[-1, :]
spitad = mpitad[-1, :]

mgbm = mgbm[:-1, :]
mmng = mmng[:-1, :]
mpitad = mpitad[:-1, :]



df = pd.DataFrame()
df["gbm_sum"] = np.sum(mgbm, axis=0)
df["gbm_cnt"] = mgbm.shape[0]
df["mng_sum"] = np.sum(mmng, axis=0)
df["mng_cnt"] = mmng.shape[0]
df["pitad_sum"] = np.sum(mpitad, axis=0)
df["pitad_cnt"] = mpitad.shape[0]
df["pgbm"] = df.gbm_sum / df.gbm_cnt
df["pmng"] = df.mng_sum / df.mng_cnt
df["ppitad"] = df.pitad_sum / df.pitad_cnt

df["smp"] = sgbm.astype(np.int32)
df["smp"] = smng.astype(np.int32)
df["smp"] = spitad.astype(np.int32)

np.linalg.norm(df.smp - df.pgbm)
np.linalg.norm(df.smp - df.pmng)
np.linalg.norm(df.smp - df.ppitad)

np.sum(df.smp == np.round(df.pgbm)) / len(df.smp)
np.sum(df.smp == np.round(df.pmng)) / len(df.smp)
np.sum(df.smp == np.round(df.ppitad)) / len(df.smp)


np.sum((df.pgbm == df.pmng) & (df.pmng == df.ppitad))
