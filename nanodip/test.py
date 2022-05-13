import bisect
from minknow_api.tools import protocols
import numpy as np
import logging
import os
import argparse
import grpc
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
    ANNOTATIONS_ABBREVIATIONS_BASEL,
    ANNOTATIONS_ABBREVIATIONS_TCGA,
    BARCODE_NAMES,
    BETA_VALUES,
    BROWSER_FAVICON,
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
    RELEVANT_GENES,
    RESULT_ENDING,
    SAMTOOLS,
    THIS_HOST,
    UMAP_PLOT_TOP_MATCHES,
)
from utils import (
    date_time_string_now,
    files_by_ending,
    discrete_colors,
    composite_path,
)
from data import (
    get_sample_methylation,
    Sample,
    Reference,
    Genome,
    get_reference_methylation,
)
from plots import (
    CNVData,
    UMAPData,
    pie_chart,
)
from webui import (
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
from api import (
    connection_from_device_id,
    parse_args,
    is_position_selected,
    predominant_barcode,
    methylation_caller,
)
import config, data, plots, nanodip, utils, api

# Define logger
logger = logging.getLogger(__name__)


sample_name = "B2021_48459_20211112_BC10"
sample_name = "B2021_48700_20211112_BC11"
sample_name = "test28"
reference_name = "GSE90496_IfP01"
reference_name = "AllIDATv2_20210804"

sample = Sample(sample_name)
reference = Reference(reference_name)
# genome = Genome()

# data.make_binary_reference_data()
cnv = CNVData(sample_name)

umapp = UMAPData(sample_name, reference_name)
umapp.read_from_disk()

umapp.sample = Sample(sample_name)
umapp.reference = Reference(reference_name)
umapp.sample.set_cpg_overlap(reference)

# sentrix_id = "201869680197_R07C01"
sentrix_id = "9968646165_R01C02"
reference_umap = "UMAP_all_bVals_top_25000.xlsx"
reference_id = "GSE90496_IfP01"
device_id = "MN26636"


# download_epidip_data(sentrix_id, reference_umap)
# umap_data = UMAPData(sentrix_id, reference_id)
# umap_data.read_precalculated_umap_matrix(reference_umap)
# umap_data.draw_pie_chart()
# umap_data.draw_scatter_plots()

import math

obs = 17
g_len = 189060
G_len = CNVData.genome.length

p_0 = g_len / G_len
N = 11161


def binomial_ci_wilson(hits, trials):
    """Return 0.99 conficence intervall of binomial distribution
    using inexact Wilson method.
    """
    p_hat = hits / trials
    return [
        p_hat - 2.58 * math.sqrt(p_hat * (1 - p_hat) / trials),
        p_hat + 2.58 * math.sqrt(p_hat * (1 - p_hat) / trials),
    ]


def _extreme_cn(observed, gene_length, sample):
    """Returns true iff expected copy number is outside of
    binomial confidence intervall of observed values.
    """
    p_0 = gene_length / CNVData.genome.length
    if not sample.reads:
        sample.set_reads()
    trials = len(sample.reads)
    p_low, p_up = binomial_ci_wilson(observed, trials)
    return p_low > p_0 or p_0 > p_up


from scipy.stats import binom


def extreme_cn(observed, gene_length, sample):
    """Returns true if and only if observed copy number is outside of
    0.99 binomial confidence interval.
    """
    p_hit = gene_length / len(CNVData.genome)
    if not sample.reads:
        sample.set_reads()
    trials = len(sample.reads)
    lower = binom.ppf(0.001, trials, p_hit)
    upper = binom.ppf(0.999, trials, p_hit)
    # lower = trials*p_hit - 3.29*math.sqrt(trials*p_hit*(1 - p_hit))
    # upper = trials*p_hit + 3.29*math.sqrt(trials*p_hit*(1 - p_hit))
    return lower > observed or observed > upper


def binom_p_value(observed, gene_length, sample):
    p_hit = gene_length / len(CNVData.genome)
    if not sample.reads:
        sample.set_reads()
    trials = len(sample.reads)
    cdf = binom.cdf(observed, trials, p_hit)
    return cdf if cdf < 0.5 else (1 - cdf)


cnv.read_from_disk()

genes = cnv.genes

# genes["extreme_cn"] = genes.apply(
# lambda z: extreme_cn(z["cn_obs"], z["len"], sample),
# axis=1,
# )

# genes["cn_p_value"] = genes.apply(
# lambda z: binom_p_value(z["cn_obs"], z["len"], sample),
# axis=1,
# )
# genes[genes.extreme_cn == True]

binom.cdf(20, 70, 0.3083573487)
binom.cdf(20, 70, 0.3083573487)
p_low, p_up = binomial_ci_wilson(obs, N)
b = extreme_cn(obs, g_len, sample)

p = 0.7
N = 100
binom.ppf(0.025, N, p)
N * p - 1.96 * math.sqrt(N * p * (1 - p))
