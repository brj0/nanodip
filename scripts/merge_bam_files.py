import os
import re
import shutil
import sys
import time

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

sys.path.insert(0, "/applications/nanodip")

from nanodip.config import (
    ENDING,
    NANODIP_OUTPUT,
)

from nanodip.utils import (
    files_by_ending,
)

from nanodip.data import (
    Genome,
    Reference,
    Sample,
)

ANNOTATION_ID = "1qmis4MSoE0eAMMwG6xZbDCrs-F1jXECZvc4wKWLR0KY"
CHROM = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
    "chr1_gl000191_random",
    "chr1_gl000192_random",
    "chr4_ctg9_hap1",
    "chr4_gl000193_random",
    "chr4_gl000194_random",
    "chr6_cox_hap2",
    "chr6_dbb_hap3",
    "chr6_mann_hap4",
    "chr6_mcf_hap5",
    "chr6_qbl_hap6",
    "chr6_ssto_hap7",
    "chr7_gl000195_random",
    "chr9_gl000198_random",
    "chr9_gl000199_random",
    "chr11_gl000202_random",
    "chr17_ctg5_hap1",
    "chr17_gl000203_random",
    "chr17_gl000204_random",
    "chr17_gl000205_random",
    "chr17_gl000206_random",
    "chr19_gl000208_random",
    "chr19_gl000209_random",
    "chrUn_gl000211",
    "chrUn_gl000212",
    "chrUn_gl000213",
    "chrUn_gl000214",
    "chrUn_gl000216",
    "chrUn_gl000217",
    "chrUn_gl000218",
    "chrUn_gl000219",
    "chrUn_gl000220",
    "chrUn_gl000221",
    "chrUn_gl000222",
    "chrUn_gl000224",
    "chrUn_gl000225",
    "chrUn_gl000226",
    "chrUn_gl000227",
    "chrUn_gl000228",
    "chrUn_gl000229",
    "chrUn_gl000230",
    "chrUn_gl000231",
    "chrUn_gl000232",
    "chrUn_gl000233",
    "chrUn_gl000234",
    "chrUn_gl000235",
    "chrUn_gl000237",
    "chrUn_gl000238",
    "chrUn_gl000239",
    "chrUn_gl000240",
    "chrUn_gl000241",
    "chrUn_gl000243",
    "chrUn_gl000246",
]
DIR = "/data/nanodip_output/merged_class/%s"
SHEET_NAME = "Sample%20list"
ANNOTATION_URL = (
    "https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s"
)


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


def annotated_cases():
    """Returns data frame containing all cases with id, b-number and
    annotation.
    """
    # Dirs containing nanodip output
    nanodip_files = [
        x for x in os.listdir(NANODIP_OUTPUT) if x.startswith("B")
    ]

    # Annotation spreadsheed
    annotation = pd.read_csv(
        ANNOTATION_URL % (ANNOTATION_ID, SHEET_NAME),
        usecols=["Sample_Name", "Methylation\nClass"],
    )
    annotation.columns = ["id_", "meth_grp"]
    annotation.id_ = annotation.apply(
        lambda x: parse_case_nr(x["id_"]), axis=1
    )
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
    return cases_df


def merge_cases(cases, fname):
    """Merges all bam files that are within the filetree of the
    files in {cases.id_}.
    """
    all_files = []

    for c in cases.id_.to_list():
        all_files.extend(files_by_ending(NANODIP_OUTPUT, c, "bam"))

    chunk_size = 50
    out = DIR % (fname + ".bam")
    tmp = DIR % "temp.csv"

    for index in tqdm(range(0, len(all_files), chunk_size)):
        files_ = all_files[index : (index + chunk_size)]
        if index > 0:
            files_ = [tmp] + files_
        pysam.merge("-fo", out, *files_)
        shutil.copyfile(out, tmp)

    os.remove(tmp)
    pysam.index(out)


def methylation_atlas(cases, fname):
    """Merges all freq_tsv ending files that are within the filetree
    of the files in {cases.id_}.
    """
    all_files = []

    for c in cases.id_.to_list():
        all_files.extend(
            files_by_ending(NANODIP_OUTPUT, c, ENDING["freq_tsv"])
        )

    methyl_df_list = []

    for file_ in tqdm(all_files):
        next_df = pd.read_csv(file_, sep="\t")
        methyl_df_list.append(next_df)

    methyl_df = pd.concat(methyl_df_list)[
        ["chromosome", "start", "methylated_frequency"]
    ]
    methyl_df["nreads"] = None
    methyl_df["methylated"] = methyl_df.methylated_frequency
    methyl_df = (
        methyl_df.groupby(["chromosome", "start"])
        .agg({"nreads": "size", "methylated": "sum"})
        .reset_index()
    )
    methyl_df["methylated_frequency"] = methyl_df.methylated / methyl_df.nreads
    chrom_to_idx = {c: (i + 1) for i, c in enumerate(CHROM)}
    methyl_df["chrom_nr"] = methyl_df.chromosome.apply(
        lambda x: chrom_to_idx[x]
    )
    methyl_df = methyl_df.sort_values(by=["chrom_nr", "start"])
    methyl_df = methyl_df.reset_index(drop=True)

    methyl_df.to_csv(DIR % (fname + ".csv"), index=False)
    methyl_df.to_pickle(DIR % (fname + ".pickle"))


class MethylAtlas:
    """Container for Methylation Atlas containing all methylation
    sites and providing search member functions.
    """

    chrom_to_idx = {c: (i + 1) for i, c in enumerate(CHROM)}

    def __init__(self, df):
        self.df = df
        self.search_idx = np.rec.fromarrays([df.chrom_nr, df.start])

    @classmethod
    def from_pickle(cls, name):
        df = pd.read_pickle(DIR % (name + ".pickle"))
        return cls(df)

    def idx(self, query_chrom, query_start):
        """Binary index search."""
        query = (MethylAtlas.chrom_to_idx[query_chrom], query_start)
        return np.searchsorted(
            self.search_idx, np.array(query, dtype=self.search_idx.dtype)
        )

    def methyl(self, left_chrom, left_start, right_chrom, right_start):
        """Returns methylation sites within an interval."""
        idx0 = self.idx(left_chrom, left_start)
        idx1 = self.idx(right_chrom, right_start)
        if (
            self.df.iloc[idx1].chromosome == right_chrom
            and self.df.iloc[idx1].start == right_start
        ):
            idx1 += 1
        return self.df[idx0:idx1]


# 17 cases
cases_df = annotated_cases()
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

# merge_cases(gbm_all[:20], "20_gbm")
# merge_cases(mng_all[:20], "20_mng")
# merge_cases(pitad_all[:20], "20_pitad")

methylation_atlas(gbm_all[:20], "20_gbm_methyl")

gbm_atlas = MethylAtlas.from_pickle("20_gbm_methyl")
gbm_atlas.methyl("chr1", 10468, "chr2", 2)
gbm_atlas.df[1:807848]

analysis_sample = gbm_all.iloc[20].id_
