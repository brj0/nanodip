import pandas as pd
import re

import sys
sys.path.insert(0, "/applications/nanodip")

from nanodip.data import Reference

reference_name = "MNG_IfP_v1"
reference_name = "AllIDATv2_20210804"
reference_name = "GSE90496_IfP01"

supervised_clfs = pd.read_csv(f"~/Documents/clf_data/{reference_name}/all_clf.csv")

URL = "https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s"
ANNOTATION_ID = "1qmis4MSoE0eAMMwG6xZbDCrs-F1jXECZvc4wKWLR0KY"
SHEET_NAME = "Sample%20list"
TRANSLATION_ID = "1hUnyA8axrJf0ppZ5oHl0nAWZJFN9obMoiPNsL3_F8_c"
TRANSLATION_NAME = "landing_zone"

annotation_raw = pd.read_csv(URL % (ANNOTATION_ID, SHEET_NAME))
translation = pd.read_csv(URL % (TRANSLATION_ID, TRANSLATION_NAME))
name_to_translation = {
    u: t for u, t in zip(translation.unique_ms_01, translation.translation)
}

annotation = pd.DataFrame(
    annotation_raw[["Sample_Name", "Methylation\nClass"]]
)
annotation.columns = ["case", "meth_grp"]

def parse_case_nr(name):
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


supervised_clfs.case = supervised_clfs.apply(
    lambda x: parse_case_nr(x["case"]), axis=1
)

annotation.case = annotation.apply(lambda x: parse_case_nr(x["case"]), axis=1)
annotation = annotation.loc[~annotation.meth_grp.isnull()]
case_to_meth = {c: m for c, m in zip(annotation.case, annotation.meth_grp)}


def get_meth_class(case_nr):
    try:
        result = case_to_meth[case_nr]
    except KeyError:
        result = "NOT_FOUND"
    return result


def common_translation(name):
    try:
        return name_to_translation[name]
    except KeyError:
        print(name)
        return name

reference = Reference(reference_name)
reference.translated_meth_grp = [
    common_translation(x) for x in reference.methylation_class
]

supervised_clfs["meth_grp"] = supervised_clfs.apply(
    lambda x: get_meth_class(x["case"]), axis=1
)

#count occurences of methylation groups in reference
meth_grp_cnt = {m:0 for m in reference.methylation_class}
for x in reference.methylation_class:
    meth_grp_cnt[x] += 1

# Methylation class must be relevant in reference
supervised_clfs["relevant"] = [
    (meth_grp_cnt.get(x, 0) > 0) for x in supervised_clfs.meth_grp
]

meth_sup_list = [f"meth_grp{d}" for d in range(1, 11)]
meth_all_list = meth_sup_list + ["meth_grp"]
prob_list = [f"prob{d}" for d in range(1, 11)]

supervised_clfs[meth_all_list] = supervised_clfs[meth_all_list].applymap(
    common_translation
)

rotated_columns = (
    ["case", "clf", "acc", "time", "relevant", "meth_grp"] + meth_sup_list + prob_list
)

# Rotate columns
supervised_clfs = supervised_clfs[rotated_columns]

# Remove duplicates
supervised_clfs = supervised_clfs.loc[
    ~supervised_clfs[["case", "clf"]].duplicated()
]

# Remove non valid cases
supervised_clfs = supervised_clfs.loc[
    ~supervised_clfs.meth_grp.isin(["-", "0_DIAGNOSTICS", "NOT_FOUND"])
]

# Sort
supervised_clfs = supervised_clfs.sort_values(by=["case", "clf"])

supervised_clfs.to_csv(f"~/Documents/clf_data/{reference_name}/all_clf_vs_annotation.csv")


# Extract cases relevant for reference data
supervised_clfs = supervised_clfs.loc[supervised_clfs.relevant]


clfs = ["RandomForestClassifier", "KNeighborsClassifier", "MLPClassifier"]
print()
print(reference_name)
print("Correctly classified cases (exact match)")
print("----------------------------------------")
for clf in clfs:
    df_clf = supervised_clfs.loc[supervised_clfs.clf == clf]
    nr_rf = len(df_clf)
    nr_rf_correct = len(df_clf.loc[df_clf.meth_grp == df_clf.meth_grp1])
    print(
        clf,
        f"{nr_rf_correct}/{nr_rf}",
        round(nr_rf_correct / nr_rf * 100, 2),
        "%",
    )

print()
print("Correctly classified cases (within top 10)")
print("------------------------------------------")
for clf in clfs:
    df_clf = supervised_clfs.loc[supervised_clfs.clf == clf]
    nr_rf = len(df_clf)
    nr_rf_correct = len(
        df_clf.loc[
            (df_clf.meth_grp == df_clf.meth_grp1)
            | (df_clf.meth_grp == df_clf.meth_grp2)
            | (df_clf.meth_grp == df_clf.meth_grp3)
            | (df_clf.meth_grp == df_clf.meth_grp4)
            | (df_clf.meth_grp == df_clf.meth_grp5)
            | (df_clf.meth_grp == df_clf.meth_grp6)
            | (df_clf.meth_grp == df_clf.meth_grp7)
            | (df_clf.meth_grp == df_clf.meth_grp8)
            | (df_clf.meth_grp == df_clf.meth_grp9)
            | (df_clf.meth_grp == df_clf.meth_grp10)
        ]
    )
    print(
        clf,
        f"{nr_rf_correct}/{nr_rf}",
        round(nr_rf_correct / nr_rf * 100, 2),
        "%",
    )

df_comb = pd.pivot(
    supervised_clfs,
    index="case",
    columns="clf",
    values=["meth_grp", "meth_grp1"],
)
df_comb = df_comb[df_comb.columns[[0, 3, 4, 5]]]
df_comb.columns = ["meth_grp", "knn", "mlp", "rf"]

print()
print("Correctly classified cases (exact match for at least 1 classifier)")
print("------------------------------------------------------------------")
nr_all = len(df_comb)
nr_all_correct = len(
    df_comb.loc[
        (df_comb.meth_grp == df_comb.knn)
        | (df_comb.meth_grp == df_comb.mlp)
        | (df_comb.meth_grp == df_comb.rf)
    ]
)
print(
    "Combined clf",
    f"{nr_all_correct}/{nr_all}",
    round(nr_all_correct / nr_all * 100, 2),
    "%",
)
print()
