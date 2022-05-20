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
    REFERENCE_METHYLATION_SHAPE,
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
    reference_methylation_from_index,
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
reference_name = "AllIDATv2_20210804"
reference_name = "GSE90496_IfP01"

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



binom.cdf(20, 70, 0.3083573487)
binom.cdf(20, 70, 0.3083573487)
p_low, p_up = binomial_ci_wilson(obs, N)
b = extreme_cn(obs, g_len, sample)

p = 0.7
N = 100
binom.ppf(0.025, N, p)
N * p - 1.96 * math.sqrt(N * p * (1 - p))

with open(REFERENCE_METHYLATION_SHAPE, "r") as f:
    num_ref, num_cpgs = [int(s) for s in f.read().splitlines()]

from sklearn.ensemble import RandomForestClassifier

sample.set_cpg_overlap(reference)
X = reference_methylation_from_index(reference.specimens_index, sample.cpg_overlap_index)
x_sample = get_sample_methylation(sample, reference)
y_raw = reference.methylation_class

mc_to_int = {}
int_to_mc = {}

for i, mc in enumerate(set(y_raw)):
    mc_to_int[mc] = i
    int_to_mc[i] = mc
y = [mc_to_int[i] for i in y_raw]

seed = 5
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.svm import SVC
X_train, X_test, y_train, y_test = train_test_split(
    X, y_raw, test_size=0.2, random_state = seed,
)

model = RandomForestClassifier(
    n_estimators=150,
    n_jobs=-1,
    random_state=seed,
)
model.fit(X_train, y_train)
model.score(X_test, y_test)


y_predict = model.predict(X_test)
# model.accuracy_score(y_test, y_predict)
y_sample_predict = model.predict([x_sample])
# y_predict =
wrong=[]
for i,j in zip(y_test,y_predict):
    if i != j:
        wrong.append((i,j))

cm = confusion_matrix(y_test, y_predict)
# fig = px.density_heatmap(cm)
# fig.write_html("/data/nanodip_reports/cm.html")
prob = model.predict_proba([x_sample])

prob_per_class = []
for p, mc in zip(prob[0], model.classes_):
    prob_per_class.append((p, mc))
prob_per_class.sort(reverse=True)

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

prob_svm = model_svm.predict_proba([x_sample])

prob_per_class_svm = []
for p, mc in zip(prob_svm[0], model_svm.classes_):
    prob_per_class_svm.append((p, mc))

def evaluate_clf(clf, x_sample, X_text, y_test):
    y_predict = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_predict)
    prob = clf.predict_proba([x_sample])[0]
    prob_per_class = []
    for p, mc in zip(prob, clf.classes_):
        prob_per_class.append((p, mc))
    prob_per_class.sort(reverse=True)
    print("Evaluation of", clf)
    print("Accuracy:", accuracy)
    for i in range(10):
        print(
            prob_per_class[i][1],
            "\t:\t",
            round(100*prob_per_class[i][0], 2),
            " %",
        )

evaluate_clf(model, x_sample, X_test, y_test)


prob_per_class_svm.sort(reverse=True)

# Finding the best hyperparameters
# params = {
    # 'C': [0.1, 1, 10, 100, 1000],
    # 'gamma': [1, 0.1, 0.01, 0.001, 0.0001],
    # 'kernel': ['rbf']
# }

# clf = GridSearchCV(
    # estimator=SVC(),
    # param_grid=params,
    # cv=5,
    # n_jobs=5,
    # verbose=1
# )

# clf.fit(X_train, y_train)
# print(clf.best_params_)

from sklearn.neighbors import KNeighborsClassifier
neigh = KNeighborsClassifier(n_neighbors=3)
neigh.fit(X_train, y_train)
print(neigh.predict([x_sample]))
prob = neigh.predict_proba([x_sample])

prob_per_class_n = []
for p, mc in zip(prob[0], neigh.classes_):
    prob_per_class_n.append((p, mc))

prob_per_class_n.sort(reverse=True)


from sklearn.neural_network import MLPClassifier
clf = MLPClassifier(solver='lbfgs', alpha=1e-5,
                    hidden_layer_sizes=(5, 2), random_state=1)

clf.fit(X_train, y_train)

print(clf.predict([x_sample]))
prob = clf.predict_proba([x_sample])
y_predict = model.predict(X_test)
accuracy_score(y_test, y_predict)
