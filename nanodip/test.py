import bisect

# import cupy
from minknow_api.tools import protocols
from multiprocessing import Pool
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
import json
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

# sys.path.insert(0, "/applications/nanodip")

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
    ILLUMINA_CPG_MAP,
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
    TMP,
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
    get_betas,
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
    single_file_methylation_caller,
    connection_from_device_id,
    remove_dirs_with_wrong_barcode,
    methylation_calling_done,
    parse_args,
    is_position_selected,
    predominant_barcode,
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

from nanodip.config import *
from nanodip.data import *
from nanodip.plots import *
from nanodip.main import *
from nanodip.utils import *
from nanodip.api import *
from nanodip.webui import *
from nanodip.classifiers import *

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
# std_all = calculate_std(reference)


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


# mgbm = get_reference_methylation(
# gbm.specimens, gbm.cpg_sites
# )
# mmng = get_reference_methylation(
# mng.specimens, mng.cpg_sites
# )
# mpitad = get_reference_methylation(
# pitad.specimens, pitad.cpg_sites
# )

# sgbm = mgbm[-1, :]
# smng = mmng[-1, :]
# spitad = mpitad[-1, :]

# mgbm = mgbm[:-1, :]
# mmng = mmng[:-1, :]
# mpitad = mpitad[:-1, :]


# df = pd.DataFrame()
# df["gbm_sum"] = np.sum(mgbm, axis=0)
# df["gbm_cnt"] = mgbm.shape[0]
# df["mng_sum"] = np.sum(mmng, axis=0)
# df["mng_cnt"] = mmng.shape[0]
# df["pitad_sum"] = np.sum(mpitad, axis=0)
# df["pitad_cnt"] = mpitad.shape[0]
# df["pgbm"] = df.gbm_sum / df.gbm_cnt
# df["pmng"] = df.mng_sum / df.mng_cnt
# df["ppitad"] = df.pitad_sum / df.pitad_cnt

# df["smp"] = sgbm.astype(np.int32)
# df["smp"] = smng.astype(np.int32)
# df["smp"] = spitad.astype(np.int32)

# np.linalg.norm(df.smp - df.pgbm)
# np.linalg.norm(df.smp - df.pmng)
# np.linalg.norm(df.smp - df.ppitad)

# np.sum(df.smp == np.round(df.pgbm)) / len(df.smp)
# np.sum(df.smp == np.round(df.pmng)) / len(df.smp)
# np.sum(df.smp == np.round(df.ppitad)) / len(df.smp)


# np.sum((df.pgbm == df.pmng) & (df.pmng == df.ppitad))


def update_num_fastq_files(sample_name):
    barcode = predominant_barcode(sample_name)
    fast5_files = [
        f
        for f in files_by_ending(DATA, sample_name, ending=".fast5")
        if barcode in f
    ]

    def from_5_to_q(fn):
        return fn.replace(".fast5", ".fastq").replace(
            "fast5_pass", "fastq_pass"
        )

    # Collect all passed fast5/fastq pairs
    fast5q_file_pairs = [
        [f, from_5_to_q(f)]
        for f in fast5_files
        if os.path.exists(from_5_to_q(f))
    ]
    num_fastq = len(fast5q_file_pairs)
    return num_fastq


def get_callable_reads_with_stats(sample_name):
    # At least 2 "passed" files need to be present.
    barcode = predominant_barcode(sample_name)
    remove_dirs_with_wrong_barcode(sample_name, barcode)
    fast5_files = [
        f
        for f in files_by_ending(DATA, sample_name, ending=".fast5")
        if barcode in f
    ]
    # Analyse in alphanumeric ordering for improved debugging.
    fast5_files.sort()

    def from_5_to_q(fn):
        return fn.replace(".fast5", ".fastq").replace(
            "fast5_pass", "fastq_pass"
        )

    # Collect all passed fast5/fastq pairs
    fast5q_file_pairs = [
        [f, from_5_to_q(f)]
        for f in fast5_files
        if os.path.exists(from_5_to_q(f))
    ]

    f5c_analysis_dir = os.path.join(NANODIP_OUTPUT, sample_name)
    if not os.path.exists(f5c_analysis_dir):
        os.mkdir(f5c_analysis_dir)

    prev_called = []
    not_called = []

    for f5, fq in fast5q_file_pairs:
        file_name = os.path.basename(f5).split(".")[0]
        analysis_dir = os.path.join(f5c_analysis_dir, file_name)
        symlink5 = os.path.join(analysis_dir, file_name + ".fast5")
        symlinkq = os.path.join(analysis_dir, file_name + ".fastq")
        if not os.path.exists(analysis_dir):
            os.mkdir(analysis_dir)
        if not os.path.exists(symlink5):
            os.symlink(f5, symlink5)
        if not os.path.exists(symlinkq):
            os.symlink(fq, symlinkq)
        if methylation_calling_done(analysis_dir):
            prev_called.append(file_name)
        else:
            not_called.append([analysis_dir, file_name])

    stats = {
        "barcode": barcode,
        "num_completed": len(prev_called),
        "num_fastq": len(fast5q_file_pairs),
        # "no_callable_left": num_fastq == num_completed,
        # "time": date_time_string_now(),
    }
    return (not_called, stats)


def methylation_caller(sample_name, analyze_one=True):
    """Searches for callable fast5/fastq files that have not yet been
    called and invokes methylation calling. Results will be added to
    the NANODIP_OUTPUT directory.

    Args:
        sample_name: Name of sample to be analyzed.
        analyse_one: If True only first fast5/fastq file found
                     will be analyzed.
    """
    not_called, stats = get_callable_reads_with_stats(sample_name)
    curr_called = []
    for directory, file_name in not_called:
        single_file_methylation_caller(directory)
        curr_called.append(file_name)
        if analyze_one:
            break
    stats["num_completed"] += len(curr_called)
    stats["no_callable_left"] = stats["num_fastq"] == stats["num_completed"]
    stats["called"] = curr_called
    stats["time"] = date_time_string_now()
    return stats


# sample_name = "test1123b"
# methylation_caller("test1123b", max_calls=0)
# methylation_caller("test1122", max_calls=0)

# def methylation_caller(sample_name, analyze_one=True):
# """Searches for callable fast5/fastq files that have not yet been
# called and invokes methylation calling. Results will be added to
# the NANODIP_OUTPUT directory.

# Args:
# sample_name: Name of sample to be analyzed.
# analyse_one: If True only first fast5/fastq file found
# will be analyzed.
# """
# # At least 2 "passed" files need to be present.
# barcode = predominant_barcode(sample_name)
# remove_dirs_with_wrong_barcode(sample_name, barcode)
# fast5_files = [
# f for f in files_by_ending(DATA, sample_name, ending=".fast5")
# if barcode in f
# ]
# # Analyse in alphanumeric ordering for improved debugging.
# fast5_files.sort()
# def from_5_to_q(fn):
# return fn.replace(
# ".fast5", ".fastq"
# ).replace("fast5_pass", "fastq_pass")

# # Collect all passed fast5/fastq pairs
# fast5q_file_pairs = [
# [f, from_5_to_q(f)] for f in fast5_files
# if os.path.exists(from_5_to_q(f))
# ]

# f5c_analysis_dir = os.path.join(NANODIP_OUTPUT, sample_name)
# if not os.path.exists(f5c_analysis_dir):
# os.mkdir(f5c_analysis_dir)

# prev_called = []
# curr_called = []
# not_called = []

# for f5, fq in fast5q_file_pairs:
# file_name = os.path.basename(f5).split(".")[0]
# analysis_dir = os.path.join(f5c_analysis_dir, file_name)
# symlink5 = os.path.join(analysis_dir, file_name + ".fast5")
# symlinkq = os.path.join(analysis_dir, file_name + ".fastq")
# if not os.path.exists(analysis_dir):
# os.mkdir(analysis_dir)
# if not os.path.exists(symlink5):
# os.symlink(f5, symlink5)
# if not os.path.exists(symlinkq):
# os.symlink(fq, symlinkq)
# if methylation_calling_done(analysis_dir):
# prev_called.append(file_name)
# else:
# not_called.append(
# [analysis_dir, file_name]
# )
# for directory, file_name in not_called:
# single_file_methylation_caller(directory)
# curr_called.append(file_name)
# if analyze_one:
# break
# num_completed = len(prev_called) + len(curr_called)
# num_fastq = len(fast5q_file_pairs)
# no_callable_left = num_fastq == num_completed
# return {
# "barcode": barcode,
# "called": curr_called,
# "num_completed": num_completed,
# "num_fastq": num_fastq,
# "no_callable_left": no_callable_left,
# "time": date_time_string_now(),
# }


# print("0         smp =", smp)
# print("1         curr_stats =", curr_stats)
# print("2         old_stats =", old_stats)
# print("3         curr_stats =", curr_stats)
# print("4         status =", status)
# print("5         progress =", progress)
# print("6         (smp, progress[smp]) =", (smp, progress[smp]))
# if (not (smp, progress[smp]) in [
# ("test1204a", "184 / 184"),
# ("test1201d", "1 / 60"),
# ("test1201c", "0 / 0"),
# ("test1201b", "0 / 92"),
# ("test1201a", "6 / 107"),
# ("test1130b", "0 / 0"),
# ("test1130a", "16 / 17"),
# ("test1129a", "12 / 19"),
# ]) and (smp != "test1201a"):
# # import pdb; pdb.set_trace()
# print("0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# x=0
# y=1
# z=2

def _get_betas(reference):
    if gpu_enabled():
        # Get unified pool
        pool = cupy.cuda.MemoryPool(cupy.cuda.memory.malloc_managed)
        # Set unified pool as default allocator
        cupy.cuda.set_allocator(pool.malloc)
        # Release GPU memory
        pool.free_all_blocks()
    else:
        logger.info("Probably no CUDA device present.")
    cpg_index = [Reference.cpg_site_to_index[c] for c in cpgs]
    specimen_bin_files = [
        composite_path(BETA_VALUES, s, ENDING["betas_bin"])
        for s in reference_specimens
    ]
    specimens_cnt = len(reference_specimens)
    cpg_cnt = len(reference.cpg_sites)
    # Determine size of beta value array. Number of CpGs is typically fixed,
    # Number of cases is variable
    block_size = GPU_RAM_USAGE // (GPU_FLOAT_SIZE * len(reference_specimens))
    beta_values_xp = xp.empty(
        [specimens_cnt, cpg_cnt], dtype=np.float64, order="C"
    )
    # Break data into blocks along the cpg-columns, adjusted to
    # GPU RAM availability.
    for col0 in tqdm(range(0, cpg_cnt, block_size), desc="Reading binary files"):
        col1 = min(cpg_cnt, col0 + block_size)
        # Will be =block_size except for the last run.
        d_col = col1 - col0
        col_index = [x-col0 for x in cpg_index if x >= col0 and x <= col1]
        for idx, file_ in enumerate(specimen_bin_files):
            beta_values_xp[idx][:d_col] = xp.fromfile(
                file_, count=d_col, offset=col0, dtype=np.float64
            )
        # Replace nan with 0.49
        beta_values_xp = xp.nan_to_num(beta_values_xp, nan=0.49)
    # Convert cupy to numpy array
    if isinstance(beta_stds_xp, cupy.ndarray):
        beta_values = beta_stds_xp.get()
    else:
       beta_values= beta_values_xp
    # Need to release GPU memory explicitly
    del beta_stds_xp
    # Release GPU memory
    if gpu_enabled():
        pool.free_all_blocks()
    return beta_values


"""
EPIC450K komplett
hg19_cnv : grid.json
reference_annotation : acronyms, all csv files
"""


def cachify_path(file_path):
    """Returns the transformed file path within the cache directory."""
    cache_name = os.path.split(CACHE)[1]
    path_components = file_path.split("/")
    path_components.insert(1, cache_name)
    new_path = "/".join(path_components)
    return new_path


# REFERENCE_METHYLATION_SHAPE = cachify_path(REFERENCE_METHYLATION_SHAPE)







reference = Reference("AllIDATv2_20210804_HPAP_Sarc")
reference_specimens = reference.specimens
cpgs = top_variable_cpgs(reference, 50000)
betas = get_betas(reference_specimens, cpgs)


