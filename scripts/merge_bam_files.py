import os
import re
import shutil
import sys

import pandas as pd
import pysam
from tqdm import tqdm

sys.path.insert(0, "/applications/nanodip")

from nanodip.config import (
    NANODIP_OUTPUT,
)

from nanodip.utils import (
    files_by_ending,
)

URL = "https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s"
ANNOTATION_ID = "1qmis4MSoE0eAMMwG6xZbDCrs-F1jXECZvc4wKWLR0KY"
SHEET_NAME = "Sample%20list"


def parse_case_nr(name):
    """Extracts valid B-number from string name.

    Example:
        >>> parse_case_nr("B2017.26170_1a")
        'B2017_26170'
    """
    name = str(name)
    try:
        result = re.search("^[A-Z]+[\d]+.[\d]+", name).group(0)
    except AttributeError:
        return name
    result = result.replace(".", "_")
    if result[3] == "_":
        print(f"Typo found in {name}")
        result = "B20" + result[1:]
    return result


# Dirs containing nanodip output
nanodip_files = [x for x in os.listdir(NANODIP_OUTPUT) if x.startswith("B")]

# Annotation spreadsheed
annotation = pd.read_csv(
    URL % (ANNOTATION_ID, SHEET_NAME),
    usecols=["Sample_Name", "Methylation\nClass"],
)
annotation.columns = ["id_", "meth_grp"]
annotation.id_ = annotation.apply(lambda x: parse_case_nr(x["id_"]), axis=1)
annotation = annotation[~annotation.meth_grp.isnull()]
annotation = annotation[annotation.meth_grp != "0_DIAGNOSTICS"]
annotation = annotation[annotation.meth_grp != "0_TRAINING"]
annotation = annotation[annotation.meth_grp != "-"]
annotation = annotation[~annotation.id_.duplicated(keep=False)]

case_to_methgrp = {x.id_: x.meth_grp for _, x in annotation.iterrows()}

# Data frame for all nanopore files
cases_df = pd.DataFrame(columns=["id_", "b_nr", "meth_grp"])
cases_df.id_ = nanodip_files
cases_df.b_nr = cases_df.apply(lambda x: parse_case_nr(x["id_"]), axis=1)
cases_df.meth_grp = cases_df.apply(
    lambda x: case_to_methgrp.get(x["b_nr"], "-"), axis=1
)


def merge_cases(cases, fname):
    """Merges all bam files that are within the filetree of the
    files in {cases.id_}.
    """
    all_files = []

    for c in cases.id_.to_list():
        all_files.extend(files_by_ending(NANODIP_OUTPUT, c, "bam"))

    chunk_size = 50
    out = "/data/nanodip_output/merged_class/%s.bam" % fname
    tmp = "/data/nanodip_output/merged_class/tmp.bam"

    for index in tqdm(range(0, len(all_files), chunk_size)):
        files_ = all_files[index : (index + chunk_size)]
        if index > 0:
            files_ = [tmp] + files_
        pysam.merge("-fo", out, *files_)
        shutil.copyfile(out, tmp)

    os.remove(tmp)
    pysam.index(out)


# 17 cases
gbm_rtk_ii = cases_df[cases_df.meth_grp.isin(["GBM_RTK_II"])]

# 49 cases
gbm_all = cases_df[
    cases_df.meth_grp.isin(
        [
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
    )
]

# 17 cases
mng_ben_1 = cases_df[cases_df.meth_grp.isin(["MNG_BEN-1"])]

# 16 cases
mng_ben_2 = cases_df[cases_df.meth_grp.isin(["MNG_BEN-2"])]

# 58 cases
mng_all = cases_df[
    cases_df.meth_grp.isin(
        [
            "MNG_BEN-3",
            "MNG_MAL",
            "MNG_BEN-1",
            "MNG",
            "MNG_BEN-2",
            "MNG_INT-B",
            "MNG_CC",
            "MNG_INT-A",
        ]
    )
]

# 23 cases
pitad_all = cases_df[
    cases_df.meth_grp.isin(
        [
            "PITAD",
            "PITAD_FSH_LH",
            "PITAD_ACTH",
            "PITAD_STH_SPA",
            "PITAD_STH_DNS_A",
            "PITAD_TSH",
            "PITAD_STH_DNS_B",
            "PITAD_PRL",
        ]
    )
]

merge_cases(gbm_all[:20], "20_gbm")
merge_cases(mng_all[:20], "20_mng")
merge_cases(pitad_all[:20], "20_pitad")
