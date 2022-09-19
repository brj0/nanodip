import os
from tqdm import tqdm
import time

from config import (
    NANODIP_OUTPUT,
)
from classifiers import (
    fit_and_evaluate_classifiers,
)

"""
Runs supervised classifiers on all B-Numbers in /data/nanodip_output/
"""

diagnostic_b_cases = [
    x for x in os.listdir(NANODIP_OUTPUT) if x.startswith("B")
]

reference_name = "GSE90496_IfP01"
reference_name = "AllIDATv2_20210804"
reference_name = "MNG_IfP_v1"

for sample_name in tqdm(diagnostic_b_cases, desc="Calculating classifiers"):
    fit_and_evaluate_classifiers(sample_name, reference_name)
