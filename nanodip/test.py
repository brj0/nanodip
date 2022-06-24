import bisect
from minknow_api.tools import protocols
import numpy as np
import logging
import os
import argparse
import grpc
import pandas as pd
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
    bonferroni_corrected_ci,
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
from classifiers import (
    fit_and_evaluate_classifiers,
    training_test_data,
    evaluate_clf,
)

import config, data, plots, nanodip, utils, api, webui, classifiers

# Define logger
logger = logging.getLogger(__name__)


# sample_name = "B2021_48459_20211112_BC10"
# sample_name = "B2021_48700_20211112_BC11"
sample_name = "test28"
# reference_name = "AllIDATv2_20210804"
reference_name = "GSE90496_IfP01"

sample = Sample(sample_name)
reference = Reference(reference_name)
genome = Genome()

cnv = CNVData(sample_name)
cnv.read_from_disk()

# umapp = UMAPData(sample_name, reference_name)
# umapp.read_from_disk()
# umapp.sample = Sample(sample_name)
# umapp.reference = Reference(reference_name)
# umapp.sample.set_cpg_overlap(reference)

from sklearn.naive_bayes import MultinomialNB, BernoulliNB, CategoricalNB
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score

from io import StringIO
str_buffer = StringIO()


# nb_clf = SVC(kernel="linear", probability=True)
# print("Start training nb")
# nb_clf.fit(X_train, y_train)
# evaluate_clf(nb_clf, x_sample, X_test, y_test)




# cm = confusion_matrix(y_test, y_predict)
# fig = px.density_heatmap(cm)
# fig.write_html("/data/nanodip_reports/cm.html")


# params = {
    # 'n_neighbors': [3, 5, 10, 20, 50, 100],
    # 'weights': ['uniform', 'distance'],
    # 'p': [1,2,5,10],
    # 'alpha': 10.0 ** -np.arange(1, 7),
# }

# gs = GridSearchCV(
    # nn_clf,
    # knn_clf,
    # param_grid=params,
    # scoring='accuracy',
    # cv=5,
    # n_jobs=5,
    # verbose=1,
# )
# gs.fit(X_train, y_train)
# gs.best_params_

sample_name = "17_6"
reference_name = "AllIDATv2_20210804"

from sklearn import preprocessing

sample = Sample(sample_name)
reference = Reference(reference_name)
# Define training/test/sample data.
X_train, X_test, y_train, y_test = training_test_data(sample, reference)
x_sample = get_sample_methylation(sample, reference)

clf = MLPClassifier(
    verbose=True,
)
svm_linear_clf = SVC(
    kernel="linear",
    probability=True,
    verbose=True,
)
from sklearn.svm import LinearSVC

start = time.time()
clf.fit(X_train, y_train)
evaluation = evaluate_clf(clf, x_sample, X_test, y_test)
passed_time = time.time() - start

