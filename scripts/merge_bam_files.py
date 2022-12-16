"""
This script can merge bam files to a methylation atlas and compare with
sample.
"""
from operator import xor
import datetime
import os
import re
import random
import shutil
import pickle

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from nanodip.data import (
    Genome,
)

from nanodip.config import (
    ENDING,
    NANODIP_OUTPUT,
)

from nanodip.utils import (
    files_by_ending,
)

# For debugging
import sys
import IPython
import IPython.core.ultratb
sys.excepthook = IPython.core.ultratb.ColorTB()
import time
pdp = lambda x: print(x.to_string())

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
    "chr6_apd_hap1",
    "chr6_cox_hap2",
    "chr6_dbb_hap3",
    "chr6_mann_hap4",
    "chr6_mcf_hap5",
    "chr6_qbl_hap6",
    "chr6_ssto_hap7",
    "chr7_gl000195_random",
    "chr8_gl000196_random",
    "chr8_gl000197_random",
    "chr9_gl000198_random",
    "chr9_gl000199_random",
    "chr9_gl000200_random",
    "chr9_gl000201_random",
    "chr8_gl000196_random",
    "chr11_gl000202_random",
    "chr17_ctg5_hap1",
    "chr17_gl000203_random",
    "chr17_gl000204_random",
    "chr17_gl000205_random",
    "chr17_gl000206_random",
    "chr18_gl000207_random",
    "chr19_gl000208_random",
    "chr19_gl000209_random",
    "chr21_gl000210_random",
    "chr18_gl000207_random",
    "chrUn_gl000211",
    "chrUn_gl000212",
    "chrUn_gl000213",
    "chrUn_gl000214",
    "chrUn_gl000215",
    "chrUn_gl000216",
    "chrUn_gl000217",
    "chrUn_gl000218",
    "chrUn_gl000219",
    "chrUn_gl000220",
    "chrUn_gl000221",
    "chrUn_gl000222",
    "chrUn_gl000223",
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
    "chrUn_gl000236",
    "chrUn_gl000237",
    "chrUn_gl000238",
    "chrUn_gl000239",
    "chrUn_gl000240",
    "chrUn_gl000241",
    "chrUn_gl000242",
    "chrUn_gl000243",
    "chrUn_gl000244",
    "chrUn_gl000245",
    "chrUn_gl000246",
    "chrUn_gl000247",
    "chrUn_gl000248",
    "chrUn_gl000249",
]
DIR = os.path.join(
    "/data/nanodip_output/merged_class",
    datetime.datetime.now().strftime("%Y-%m-%d-%H%M"),
)
LOG_FILE = os.path.join(DIR, "methyl_atlas.log")
NUM_LOOPS = 20
NUM_LOOPS = 100  # TODO del
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

    case_to_methgrp = {x.id_: x.meth_grp for x in annotation.itertuples()}

    # Data frame for all nanopore files
    cases_df = pd.DataFrame(columns=["id_", "b_nr", "meth_grp"])
    cases_df.id_ = nanodip_files
    cases_df.b_nr = cases_df.apply(lambda x: parse_case_nr(x["id_"]), axis=1)
    cases_df.meth_grp = cases_df.apply(
        lambda x: case_to_methgrp.get(x["b_nr"], "-"), axis=1
    )
    return cases_df


def merge_bam_files(cases, fname):
    """Merges all bam files that are within the filetree of the
    files in {cases}.
    """
    bam_files = []

    for c in cases:
        bam_files.extend(files_by_ending(NANODIP_OUTPUT, c, "bam"))

    chunk_size = 50
    out = os.path.join(DIR, fname + ".bam")
    tmp = os.path.join(DIR, "tmp.bam")

    for index in tqdm(range(0, len(bam_files), chunk_size), desc="merge bam"):
        files_ = bam_files[index : (index + chunk_size)]
        if index > 0:
            files_ = [tmp] + files_
        pysam.merge("-fo", out, *files_)
        shutil.copyfile(out, tmp)

    os.remove(tmp)
    pysam.index(out)


def get_all_tsv_files(cases):
    """Returns a list containing all freq_tsv files in the filetree of
    {cases}.
    """
    files = []
    for c in cases:
        files.extend(files_by_ending(NANODIP_OUTPUT, c, ENDING["freq_tsv"]))
    return files


def merge_methylation(fname, cases=None, files=None):
    """Merges all freq_tsv ending files that are either within the
    filetree of the files in {cases} or within {files}.
    """
    if not xor(cases is None, files is None):
        raise ValueError("Provide either 'cases' or 'files'")
    all_files = files if cases is None else get_all_tsv_files(cases)
    methyl_df_list = []
    for file_ in tqdm(all_files, desc="merge tsv"):
        next_df = pd.read_csv(file_, sep="\t")
        methyl_df_list.append(next_df)
    methyl_df = pd.concat(methyl_df_list)[
        ["chromosome", "start", "methylated_frequency"]
    ]
    methyl_df["read_cnt"] = None
    methyl_df["methyl_cnt"] = methyl_df.methylated_frequency
    methyl_df = methyl_df[methyl_df.methyl_cnt.isin([0.0, 1.0])]
    methyl_df.methyl_cnt = methyl_df.methyl_cnt.astype(int)
    methyl_df = (
        methyl_df.groupby(["chromosome", "start"])
        .agg({"read_cnt": "size", "methyl_cnt": "sum"})
        .reset_index()
    )
    methyl_df["methyl_freq"] = methyl_df.methyl_cnt / methyl_df.read_cnt
    methyl_df["methylated"] = round(methyl_df.methyl_freq)
    chrom_to_idx = {c: (i + 1) for i, c in enumerate(CHROM)}
    methyl_df["chrom_nr"] = methyl_df.chromosome.apply(
        lambda x: chrom_to_idx[x]
    )
    methyl_df = methyl_df.sort_values(by=["chrom_nr", "start"])
    methyl_df = methyl_df.reset_index(drop=True)
    methyl_df.to_csv(os.path.join(DIR, fname + ".csv"), index=False)
    methyl_df.to_pickle(os.path.join(DIR, fname + ".pickle"))
    return methyl_df


class MethylClass:
    """Container for merged methylation class object containing all
    methylation sites and providing search member functions.
    """

    chrom_to_idx = {c: (i + 1) for i, c in enumerate(CHROM)}

    def __init__(self, df, name):
        self.df = df
        self.name = name
        self.search_idx = np.rec.fromarrays([df.chrom_nr, df.start])

    @classmethod
    def from_disk(cls, name):
        df = pd.read_pickle(os.path.join(DIR, name + ".pickle"))
        return cls(df, name)

    def idx(self, query_chrom, query_start):
        """Binary index search."""
        query = (MethylClass.chrom_to_idx[query_chrom], query_start)
        return np.searchsorted(
            self.search_idx, np.array(query, dtype=self.search_idx.dtype)
        )

    def methyl(self, chrom, left_start, right_start):
        """Returns methylation sites within an interval."""
        idx0 = self.idx(chrom, left_start)
        idx1 = self.idx(chrom, right_start)
        if (
            self.df.iloc[idx1].chromosome == chrom
            and self.df.iloc[idx1].start == right_start
        ):
            idx1 += 1
        return self.df[idx0:idx1]

    def __str__(self):
        """Prints overview of object for debugging purposes."""
        lines = [
            f"MethylClass object:",
            f"name: '{self.name}'",
            f"df: '{self.df}'",
            f"search_idx:\n{self.search_idx}",
        ]
        return "\n".join(lines)


def cpg_offsets(seq, offset, methyl_status):
    """Transforms DNA sequance coordinate to offset of all CpG
    sites.
    Example:
        >>> cpg_offsets('TGATCCGATCGATCG', 14850238, 1.0)
        [[14850238, 1.0], [14850242, 1.0], [14850246, 1.0]]
    """
    start0 = next(re.finditer("CG", seq)).start()
    offsets = [
        [cg.start() - start0 + offset, methyl_status]
        for cg in re.finditer("CG", seq)
    ]
    return offsets


def get_sample_reads(sample_name):
    """Return data from containing all sample reads with methylation
    status.
    """
    files_ = files_by_ending(NANODIP_OUTPUT, sample_name, ENDING["result_tsv"])
    read_list = []
    for file_ in files_:
        data = pd.read_csv(
            file_,
            delimiter="\t",
        )
        read_list.append(data)
    f5c_out = pd.concat(read_list)
    f5c_out.loc[f5c_out.log_lik_ratio >= 2.5, "methyl_status"] = 1
    f5c_out.loc[f5c_out.log_lik_ratio <= -2.5, "methyl_status"] = 0
    f5c_out = f5c_out[~f5c_out.methyl_status.isnull()]
    f5c_out["cpg_start_methyl"] = [
        cpg_offsets(x, y, z)
        for x, y, z in zip(
            f5c_out.sequence, f5c_out.start, f5c_out.methyl_status
        )
    ]
    reads = (
        f5c_out.groupby(["read_name", "chromosome"])["cpg_start_methyl"]
        .apply(sum)
        .reset_index()
        .sort_values(["chromosome"])
    ).reset_index(drop=True)
    return reads


def _atlas_similarity(sample_reads, atlas):
    """Same as atlas_similarity but loops over reads and is slower for
    large data sets.
    """
    methyl_matches = {mc: [] for mc in atlas}
    methyl_overlap = {mc: [] for mc in atlas}
    for read in tqdm(
        sample_reads.itertuples(),
        total=sample_reads.shape[0],
        desc="Compare reads with reference",
    ):
        read_methyl = read.cpg_start_methyl
        start = read_methyl[0][0]
        end = read_methyl[-1][0]
        for mclass in atlas:
            ref = mclass.methyl(read.chromosome, start, end)
            overlap = pd.merge(
                pd.DataFrame(read_methyl, columns=["start", "methylated_smp"]),
                ref[["start", "methylated"]],
                on="start",
            )
            matches = len(
                overlap[overlap.methylated_smp == overlap.methylated]
            )
            methyl_matches[mclass].append(matches)
            methyl_overlap[mclass].append(len(overlap))
    return methyl_matches, methyl_overlap


def atlas_similarity(sample_reads, atlas):
    """Returns overlapping and identical CpG methylation of sample with
    each methylation class.
    """
    sample_cpgs = sample_reads.explode("cpg_start_methyl", ignore_index=True)
    sample_cpgs[["start", "methylated_smp"]] = pd.DataFrame(
        sample_cpgs.cpg_start_methyl.to_list()
    )
    sample_cpgs = sample_cpgs.drop("cpg_start_methyl", axis=1)

    chrom_to_idx = {c: (i + 1) for i, c in enumerate(CHROM)}
    methyl_overlap = {}
    methyl_matches = {}
    cpgs = {}

    for mclass in atlas:
        cpgs_df = sample_cpgs.merge(
            mclass.df[["chromosome", "start", "methylated"]],
            how="left",
            on=["chromosome", "start"],
        )
        cpgs_df["chrom_nr"] = cpgs_df.chromosome.apply(
            lambda x: chrom_to_idx[x]
        )
        cpgs_df = cpgs_df.sort_values(by=["chrom_nr", "start"])
        cpgs_df = cpgs_df.reset_index(drop=True)
        cpgs_df["hit"] = cpgs_df.methylated_smp == cpgs_df.methylated
        cpgs_df["in_atlas"] = False
        cpgs_df.loc[~cpgs_df.methylated.isna(), "in_atlas"] = True
        summary = cpgs_df.groupby(["read_name", "chromosome"]).agg(
            hit_cnt=pd.NamedAgg(column="hit", aggfunc="sum"),
            overlap_cnt=pd.NamedAgg(column="in_atlas", aggfunc="sum"),
        )
        methyl_overlap[mclass] = summary.overlap_cnt.to_list()
        methyl_matches[mclass] = summary.hit_cnt.to_list()
        cpgs[mclass] = cpgs_df
    return methyl_matches, methyl_overlap, cpgs


def merged_atlas_with_sample(atlas, sample_reads):
    """Merge and return all methylation classes to a single data frame."""
    merged = atlas[0].df
    chrom_to_idx = {c: (i + 1) for i, c in enumerate(CHROM)}
    for mclass in tqdm(atlas[1:], desc="merge methylation classes"):
        merged = merged.merge(
            mclass.df,
            how="outer",
            on=["chromosome", "start"],
            suffixes=[None, "_" + mclass.name],
        )
    suffix0 = "_" + atlas[0].name
    merged.rename(
        columns={
            "read_cnt": "read_cnt" + suffix0,
            "methyl_cnt": "methyl_cnt" + suffix0,
            "methyl_freq": "methyl_freq" + suffix0,
            "methylated": "methylated" + suffix0,
            "chrom_nr": "chrom_nr" + suffix0,
        },
        inplace=True,
    )
    sample_cpgs = sample_reads.explode("cpg_start_methyl", ignore_index=True)
    sample_cpgs[["start", "methylated_smp"]] = pd.DataFrame(
        sample_cpgs.cpg_start_methyl.to_list()
    )
    sample_cpgs = sample_cpgs.drop("cpg_start_methyl", axis=1)
    merged = merged.merge(
        sample_cpgs,
        how="outer",
        on=["chromosome", "start"],
        suffixes=[None, "_smp"],
    )
    merged = merged[
        ~merged.methylated_gbm.isna()
        & ~merged.methylated_mng.isna()
        & ~merged.methylated_pitad.isna()
        & ~merged.methylated_smp.isna()
    ]
    merged["chrom_nr"] = merged.chromosome.apply(lambda x: chrom_to_idx[x])
    merged = merged.sort_values(by=["chrom_nr", "start"])
    merged = merged.reset_index(drop=True)
    return merged


class Run:
    """Container containing all info about a run."""

    def __init__(
        self,
        sample_name,
        methyl_matches,
        methyl_overlap,
        atlas,
        cpgs,
        reference_dict=None,
        true_class=None,
    ):
        self.sample_name = sample_name
        self.matches = methyl_matches
        self.overlap = methyl_overlap
        self.atlas = atlas
        self.cpgs = cpgs
        self.reference_dict = reference_dict
        self.true_class = true_class
        self.hit_cnt = None
        self.read_cnt = None
        self.read_cnt_net = None
        self.calculate_results()

    def calculate_results(self):
        methyl_net_matches = {
            mc: [
                2 * hit - ovlp
                for hit, ovlp in zip(self.matches[mc], self.overlap[mc])
            ]
            for mc in self.atlas
        }
        max_hit = [
            max(x) for x in zip(*(self.matches[mc] for mc in self.atlas))
        ]
        max_net_hit = [
            max(x) for x in zip(*(methyl_net_matches[mc] for mc in self.atlas))
        ]
        self.hit_cnt = {x.name: sum(self.matches[x]) for x in self.atlas}
        self.read_cnt = {
            mclass.name: sum(
                mat == max_ for mat, max_ in zip(self.matches[mclass], max_hit)
            )
            for mclass in self.atlas
        }
        self.read_cnt_net = {
            mclass.name: sum(
                mat == max_
                for mat, max_ in zip(methyl_net_matches[mclass], max_net_hit)
            )
            for mclass in self.atlas
        }

    def __str__(self):
        """Prints overview of object for debugging purposes."""
        lines = [
            f"sample_name: {self.sample_name}",
            f"true_class: {self.true_class}",
            f"hit_cnt: {self.hit_cnt}",
            f"read_cnt: {self.read_cnt}",
            f"read_cnt_net: {self.read_cnt_net}",
        ]
        return "\n".join(lines)


def make_runs(cases_all):
    min_case_per_entity = min(len(x) for x in cases_all.values())
    for iter_ in range(NUM_LOOPS):
        random_cases = {
            entity: random.sample(cases, k=min_case_per_entity)
            for entity, cases in cases_all.items()
        }
        random_samples = {
            entity: random_cases[entity][0] for entity in random_cases
        }
        random_references = {
            entity: random_cases[entity][1:] for entity in random_cases
        }

        sizes = {}
        bam_files = {}

        for entity in cases_all:
            merge_methylation(entity, random_references[entity])
            sizes[entity] = 0
            bam_files[entity] = []
            for c in random_references[entity]:
                bam_files[entity].extend(
                    files_by_ending(NANODIP_OUTPUT, c, "bam")
                )
            for f in bam_files[entity]:
                sizes[entity] += os.path.getsize(f)
            print("Size of", entity, "bam files:", sizes[entity], "bytes")

        atlas = []
        for entity in cases_all:
            atlas.append(MethylClass.from_disk(entity))

        # Restrict atlas to chr1-chr22
        for mclass in atlas:
            mclass.df = mclass.df[mclass.df.chromosome.isin(CHROM[:22])]

        # Restrict atlas to overlapping cpgs
        mclass0 = atlas[0]
        for mclass in atlas[1:]:
            mclass0.df = mclass0.df.merge(
                mclass.df[["chromosome", "start"]],
                on=["chromosome", "start"],
                suffixes=[None, "_y"],
            )
        for mclass in atlas[1:]:
            mclass.df = mclass.df.merge(
                mclass0.df[["chromosome", "start"]],
                on=["chromosome", "start"],
                suffixes=[None, "_y"],
            )

        # # Make all df's in atlas the same size
        # min_cpg_cnt = min(mclass.df.shape[0] for mclass in atlas)
        # for mclass in atlas:
        # rand_idx = random.sample(range(mclass.df.shape[0]), k=min_cpg_cnt)
        # rand_idx.sort()
        # mclass.df = mclass.df.loc[rand_idx]

        # Print atlases
        for mclass in atlas:
            print(mclass)

        for entity in cases_all:
            sample_name = random_samples[entity]
            sample_reads = get_sample_reads(sample_name)
            sample_reads["true_class"] = entity
            methyl_matches, methyl_overlap, cpgs = atlas_similarity(
                sample_reads,
                atlas,
            )
            run = Run(
                sample_name,
                methyl_matches,
                methyl_overlap,
                atlas,
                cpgs,
                random_references,
                entity,
            )
            with open(LOG_FILE, "a") as f:
                f.write(str(run) + "\n\n")
            with open(
                os.path.join(DIR, f"{entity}_{iter_}_run.pickle"), "wb"
            ) as f:
                pickle.dump(run, f)


def evaluate_runs():
    run_list = []
    pickle_files = files_by_ending(DIR, "", "run.pickle")
    for file_ in tqdm(pickle_files, desc="reading pickle files"):
        with open(file_, "rb") as f:
            run = pickle.load(f)
        run_list.append(
            {
                "sample_name": run.sample_name,
                "true_class": run.true_class,
                "hit_cnt": run.hit_cnt,
                "read_cnt": run.read_cnt,
                "read_cnt_net": run.read_cnt_net,
            }
        )
        del run

    # Evaluate classification
    for atlas_nm in ["gbm", "mng", "pitad"]:
        print("\n# " + atlas_nm)
        for method in ["hit_cnt", "read_cnt", "read_cnt_net"]:
            correct = [
                r
                for r in run_list
                if r["true_class"] == atlas_nm
                and r[method][atlas_nm] == max(r[method].values())
            ]
            print(f"{method}: {len(correct)}/{len(run_list)/3}")

    # Print number of cpgs in atlases for each chromosome
    print("# CpG distribution of atlases")
    for chrom in CHROM:
        line = [chrom, "\t"]
        for mclass in atlas:
            line.extend(
                [
                    mclass.name,
                    "=",
                    len(mclass.df[mclass.df.chromosome == chrom]),
                    " ",
                ]
            )
        print("".join([str(l) for l in line]))

    # Print number of cpgs in atlases for all chromosomes
    line = ["chr1-chr22", "\t"]
    for mclass in atlas:
        line.extend(
            [
                mclass.name,
                "=",
                len(mclass.df[mclass.df.chromosome.isin(CHROM[:22])]),
                " ",
            ]
        )

    print("# Total")
    print("".join([str(l) for l in line]))

    # Print raw data
    print("# Raw data")
    for r in run_list:
        print("sample_name", r["sample_name"])
        print("true_class", r["true_class"])
        print("hit_cnt", r["hit_cnt"])
        print("read_cnt", r["read_cnt"])
        print("read_cnt_net", r["read_cnt_net"])
        print()


def random_atlas_samples(cases_all):
    min_case_per_entity = min(len(x) for x in cases_all.values())
    random_cases = {
        entity: random.sample(cases, k=min_case_per_entity)
        for entity, cases in cases_all.items()
    }
    random_samples = {
        entity: random_cases[entity][0] for entity in random_cases
    }
    random_references = {
        entity: random_cases[entity][1:] for entity in random_cases
    }
    atlas = []
    # Merge methylation classes of atlas
    for entity in cases_all:
        mclass = merge_methylation(entity, random_references[entity])
        atlas.append(MethylClass(mclass, entity))
    # Restrict atlas to chr1-chr22
    for mclass in atlas:
        mclass.df = mclass.df[mclass.df.chromosome.isin(CHROM[:22])]
    runs = []
    for entity in cases_all:
        sample_name = random_samples[entity]
        sample_reads = get_sample_reads(sample_name)
        sample_reads["true_class"] = entity
        sample_reads["start"] = sample_reads.cpg_start_methyl.apply(
            lambda x: min(y[0] for y in x)
        )
        sample_reads["end"] = sample_reads.cpg_start_methyl.apply(
            lambda x: max(y[0] for y in x)
        )
        sample_reads["read_len"] = sample_reads.end - sample_reads.start + 1
        merged = merged_atlas_with_sample(atlas, sample_reads)
        merged["gbm_hit"] = merged.methylated_smp == merged.methylated_gbm
        merged["mng_hit"] = merged.methylated_smp == merged.methylated_mng
        merged["pitad_hit"] = merged.methylated_smp == merged.methylated_pitad
        runs.append(merged)
    return runs


# Make output dir
os.makedirs(DIR, exist_ok=True)

# 17 cases
all_cases_df = annotated_cases()
gbm_rtk_ii = all_cases_df[all_cases_df.meth_grp.isin(["GBM_RTK_II"])]

# 49 cases
gbm_all = all_cases_df[
    all_cases_df.meth_grp.isin(
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

# 33 cases
mng_ben = all_cases_df[all_cases_df.meth_grp.isin(["MNG_BEN-1", "MNG_BEN-2"])]

# 58 cases
mng_all = all_cases_df[
    all_cases_df.meth_grp.isin(
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
pitad_all = all_cases_df[
    all_cases_df.meth_grp.isin(
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

# 12 cases
pitad_fsh = all_cases_df[
    all_cases_df.meth_grp.isin(
        [
            "PITAD_FSH_LH",
        ]
    )
]

cases_all = {
    "gbm": gbm_all.id_.tolist(),
    "mng": mng_ben.id_.tolist(),
    "pitad": pitad_all.id_.tolist(),
}

# make_runs(cases_all)


atlas_samples = []
for iter_ in tqdm(range(NUM_LOOPS), desc="sampling loops"):
    rand_samp = random_atlas_samples(cases_all)
    merged_rand_samp = pd.concat(rand_samp)
    merged_rand_samp["iter"] = iter_
    atlas_samples.append(merged_rand_samp)

merged = pd.concat(atlas_samples)

summary = (
    merged.groupby(
        [
            "read_name",
            "chromosome",
            "read_len",
            "chrom_nr",
            "iter",
            "true_class",
        ]
    )
    .agg(
        overlap_cnt=pd.NamedAgg(column="methylated_smp", aggfunc="size"),
        gbm=pd.NamedAgg(column="gbm_hit", aggfunc="sum"),
        mng=pd.NamedAgg(column="mng_hit", aggfunc="sum"),
        pitad=pd.NamedAgg(column="pitad_hit", aggfunc="sum"),
    )
    .reset_index()
)
summary = summary.sort_values(by=["chrom_nr", "overlap_cnt"])

merged.to_csv(os.path.join(DIR, "merged_cpgs.csv"), index=False)

summary["correct"] = False
all_entities = ["gbm", "mng", "pitad"]
for true_class in all_entities:
    entities = set(all_entities)
    alternatives = entities - set([true_class])
    summary.loc[
        (summary.true_class == true_class)
        & (summary[true_class] > summary[alternatives].max(axis=1)),
        "correct",
    ] = True
summary.to_csv(os.path.join(DIR, "summary_reads.csv"), index=False)

correct = []
for i in range(5000):
    summary_ = summary[summary.overlap_cnt >= i]
    if len(summary_) > 0:
        prop = sum(summary_.correct) / len(summary_)
        correct.append(prop)

np.array(correct).tofile(os.path.join(DIR, "correct_reads_5000.csv"), sep=",")

# merge_bam_files(gbm_all[:20], "20_gbm")
# merge_bam_files(gbm_all[20:21], "sample_gbm_nr20")
# merge_bam_files(mng_all[:20], "20_mng")
# merge_bam_files(pitad_all[:20], "20_pitad")

# merge_methylation("20_gbm_methyl", gbm_all[:20].id_)
# merge_methylation("sample_gbm_nr20", gbm_all[20:21].id_)
# merge_methylation("20_mng_methyl", mng_all[:20].id_)
# merge_methylation("20_pitad_methyl", pitad_all[:20].id_)
# merge_methylation("test", pitad_all[:2].id_)

# gbm_atlas = MethylClass.from_disk("20_gbm_methyl")
# mng_atlas = MethylClass.from_disk("20_mng_methyl")
# pitad_atlas = MethylClass.from_disk("20_pitad_methyl")

# sample_name = gbm_all.iloc[20].id_
# sample_reads = get_sample_reads(sample_name)
# atlas = [gbm_atlas, mng_atlas, pitad_atlas]

# methyl_matches, methyl_overlap, cpgs = atlas_similarity(
# sample_reads,
# atlas,
# )

# run = Run(sample_name, methyl_matches, methyl_overlap, atlas, cpgs)
# print(run)

# m = merged[
# [
# "chromosome",
# "start",
# "end",
# "true_class",
# "read_len",
# "read_name",
# "methylated_gbm",
# "methylated_mng",
# "methylated_pitad",
# "methylated_smp",
# ]
# ]
# readnms = list(set(merged.read_name))
# r = readnms[0]
# r = "01ed94c0-93e1-45b5-9b13-63dc3cb00926"
# m[m.read_name == r]
# read_lengths = m[
# ~m[["read_name", "chromosome"]].duplicated()
# ].read_len.to_list()
