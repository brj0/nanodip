from minknow_api.tools import protocols
from plotly.io import write_json, from_json
from scipy.stats import binomtest
from sklearn import preprocessing
from sklearn.ensemble import randomforestclassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.model_selection import gridsearchcv, train_test_split
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import multinomialnb, bernoullinb, categoricalnb
from sklearn.neighbors import kneighborsclassifier
from sklearn.neural_network import mlpclassifier
from sklearn.svm import linearsvc
from sklearn.svm import svc
import argparse
import bisect
import grpc
import math
import plotly.express as px
import plotly.graph_objects as go
import pysam
import random
import re
from io import stringio
from tqdm import tqdm
import logging
import numpy as np
import openpyxl
import os
import pandas as pd
import sys
import time

from config import (
    analysis_exclusion_patterns,
    annotations,
    annotations_abbreviations_basel,
    annotations_abbreviations_tcga,
    barcode_names,
    beta_values,
    browser_favicon,
    cherrypy_host,
    cherrypy_port,
    chromosomes,
    cnv_grid,
    cnv_link,
    data,
    debug_mode,
    ending,
    umap_link,
    epidip_umap_coordinate_files,
    excluded_from_analysis,
    f5c,
    genes,
    genes_raw,
    illumina_cg_map,
    minimap2,
    nanodip_output,
    nanodip_reports,
    needed_number_of_bases,
    plotly_render_mode,
    reads_per_file,
    reference_genome_fa,
    reference_genome_mmi,
    reference_methylation_shape,
    relevant_genes,
    result_ending,
    samtools,
    this_host,
    umap_plot_top_matches,
)

from utils import (
    date_time_string_now,
    files_by_ending,
    discrete_colors,
    composite_path,
    bonferroni_corrected_ci,
)

from data import (
    get_sample_methylation,
    sample,
    reference,
    genome,
    get_reference_methylation,
    reference_methylation_from_index,
)

from plots import (
    cnvdata,
    umapdata,
    pie_chart,
)

from webui import (
    device,
    devices,
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

from classifiers import (
    fit_and_evaluate_classifiers,
    training_test_data,
    evaluate_clf,
)


import config, data, plots, nanodip, utils, api, webui, classifiers

print("import done")

# define logger
logger = logging.getlogger(__name__)

# sample_name = "b2021_48459_20211112_bc10"
# sample_name = "b2021_48700_20211112_bc11"
# reference_name = "allidatv2_20210804"

# sample_name = "test28"
# reference_name = "gse90496_ifp01"

# sample = sample(sample_name)
# reference = reference(reference_name)
# genome = genome()

# cnv = cnvdata(sample_name)
# cnv.read_from_disk()

# umapp = umapdata(sample_name, reference_name)
# umapp.read_from_disk()
# umapp.sample = sample(sample_name)
# umapp.reference = reference(reference_name)
# umapp.sample.set_cpg_overlap(reference)

# nb_clf = svc(kernel="linear", probability=true)
# print("start training nb")
# nb_clf.fit(x_train, y_train)
# evaluate_clf(nb_clf, x_sample, x_test, y_test)

# cm = confusion_matrix(y_test, y_predict)
# fig = px.density_heatmap(cm)
# fig.write_html("/data/nanodip_reports/cm.html")

# params = {
# 'n_neighbors': [3, 5, 10, 20, 50, 100],
# 'weights': ['uniform', 'distance'],
# 'p': [1,2,5,10],
# 'alpha': 10.0 ** -np.arange(1, 7),
# }

# gs = gridsearchcv(
# nn_clf,
# knn_clf,
# param_grid=params,
# scoring='accuracy',
# cv=5,
# n_jobs=5,
# verbose=1,
# )
# gs.fit(x_train, y_train)
# gs.best_params_

reference_name ="gse90496_ifp01"
sample_name = "b1992_24268_20211126_bc12"

sample = sample(sample_name)
reference = reference(reference_name)

# define training/test/sample data.
X_train, X_test, y_train, y_test = training_test_data(sample, reference)
x_sample = get_sample_methylation(sample, reference)

clf = RandomForestClassifier(
verbose=True,
)
svm_linear_clf = SVC(
kernel="linear",
probability=True,
verbose=True,
)

start = time.time()
clf.fit(X_train, y_train)
evaluation = evaluate_clf(clf, x_sample, X_test, y_test)
passed_time = time.time() - start


y_predict = clf.predict(X_test)
# Fraction of correctly classified test samples.
accuracy = accuracy_score(y_test, y_predict)
prob = clf.predict_proba([x_sample])[0]
prob_per_class = [(p, c) for p, c in zip(prob, clf.classes_)]
prob_per_class.sort(reverse=True)




clf_knn = kneighborsclassifier(
verbose=True,
)

start = time.time()
clf_knn.fit(X_train, y_train)
evaluation = evaluate_clf(clf_knn, x_sample, X_test, y_test)
passed_time = time.time() - start


y_predict = clf_knn.predict(X_test)
# Fraction of correctly classified test samples.
accuracy = accuracy_score(y_test, y_predict)
prob = clf_knn.predict_proba([x_sample])[0]
prob_per_class = [(p, c) for p, c in zip(prob, clf_knn.classes_)]
prob_per_class.sort(reverse=True)


# #EPIDIP
# # EpiDiP functionality, temporary file directory
# EPIDIP_TMP = os.path.join(DATA, "epidip_tmp")
# GPU_ENABLED = False

# # Desired RAM usage. 4 works best on Jetson AGX 32GB; adapt to GPU / RAM
# # layout, see
# # https://forums.developer.nvidia.com/t/nvmapreserveop-0x80000000-failed-22-when-running-cufft-plans/168781/14
# GPU_RAM_USAGE = 4 * 1024**3

# # Size per float in GPU RAM (tested on AGX Xavier)
# GPU_FLOAT_SIZE = 8


# if GPU_ENABLED:
    # import cupy as xp
# else:
    # import numpy as xp


# def calculate_std(reference_name):
    # """Calculate sorted standard deviations with GPU (if present) for
    # a particular reference dataset
    # """
    # pass


# if not os.path.exists(EPIDIP_TMP):
    # os.makedirs(EPIDIP_TMP)

# std_bin = composite_path(EPIDIP_TMP, reference_name, ENDING["stdarr"])
# std_sorted_csv = composite_path(EPIDIP_TMP, reference_name, ENDING["stdsort"])

# if GPU_ENABLED:
    # # Get unified pool
    # pool = cupy.cuda.MemoryPool(cupy.cuda.memory.malloc_managed)
    # # Set unified pool as default allocator
    # cupy.cuda.set_allocator(pool.malloc)
    # # Release GPU memory
    # pool.free_all_blocks()
# else:
    # logger.info("Probably no CUDA device present.")

# # Collect reference case binary file names

# reference = Reference(reference_name)
# reference.specimens
# reference.cpg_sites

# binFilePresent = [
    # composite_path(BETA_VALUES, s, ENDING["betas"])
    # for s in reference.specimens
# ]

# # Determine size of beta value array. Number of CpGs is typically fixed,
# # number of cases is variable
# block_size = round(GPU_RAM_USAGE / GPU_FLOAT_SIZE / len(reference.specimens))
# # binFilePresent = binFilePresent[:1000] # TODO del

# # Create first column in dataframe, containing Illumina CpG names.
# beta_value_df = pd.DataFrame(reference.cpg_sites, columns=["cpg_site"])

# numFiles = len(binFilePresent)

# # rows in ilmnID column
# cpg_cnt = beta_value_df.shape[0]

# # create fixed-size cupy array filled with -1
# beta_values = xp.full([numFiles, block_size], -1, dtype=float, order="C")

# beta_stds = xp.array([])
# # Break data into blocks, adjusted to GPU RAM availability.
# for i, col0 in enumerate(tqdm(range(0, cpg_cnt, block_size), desc="Calculating std")):
    # col1 = col0 + block_size
    # if col1 > cpg_cnt:
        # col1 = cpg_cnt
        # beta_values.resize((numFiles, col1 - col0))
    # for idx, file_ in enumerate(binFilePresent):
        # idx, file_ = 0, binFilePresent[0]
        # beta_values[idx] = xp.fromfile(
            # file_, count=(col1 - col0), offset=col0, dtype=float
        # )
    # # replace nan with 0.49
    # beta_values = xp.nan_to_num(beta_values, nan=0.49)
    # beta_stds = xp.append(beta_stds, xp.std(beta_values, axis=0, dtype=float))

# # standard deviations >1 are useless (typically INF values, eliminate)
# beta_stds[beta_stds > 1] = 0

# beta_stds.tofile(std_bin)
# beta_value_df["binIndex"] = range(0, cpg_cnt)
# if str(type(beta_stds)) == "<class 'numpy.ndarray'>":
    # beta_value_df["StDev"] = beta_stds
# else:
    # beta_value_df[
        # "StDev"
    # ] = beta_stds.get()  # get is required for cupy arrays only
# beta_value_df.sort_values(
    # by="StDev",
    # axis=0,
    # ascending=False,
    # inplace=True,
    # kind="quicksort",
    # na_position="last",
    # ignore_index=False,
    # key=None,
# )
# beta_value_df.to_csv(path_or_buf=std_sorted_csv, index=False)
# # need to release GPU memory explicitly

# # TODO redo
# # del beta_value_df
# # del beta_stds
# # del beta_values

# # release GPU memory
# if GPU_ENABLED:
    # pool.free_all_blocks()
