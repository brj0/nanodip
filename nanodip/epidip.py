"""
### EpiDiP Functions
"""

# start_external_modules
from urllib import request
import logging
import os
import time

from pdf2image import convert_from_path
from tqdm import tqdm
import numpy as np
import pandas as pd

# end_external_modules

# start_internal_modules
from nanodip.config import (
    BETA_VALUES,
    CNV_LINK,
    ENDING,
    EPIDIP_TMP,
    GPU_FLOAT_SIZE,
    GPU_RAM_USAGE,
    NANODIP_REPORTS,
    UMAP_LINK,
)
from nanodip.data import (
    Reference,
)
from nanodip.utils import (
    composite_path,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)

# Import cupy if available.
try:
    import cupy as xp
    import cupy
except ImportError:
    import numpy as xp


def gpu_enabled():
    """Tests if CUDA device is present."""
    try:
        cupy.cuda.Device()
        return True
    except:
        return False


def download_epidip_data(sentrix_id, reference_umap):
    """Downloads UMAP plot coordinates of reference data and CNV plot of
    sample with given Sentrix ID.
    """
    umap_coordinates_local = composite_path(
        NANODIP_REPORTS,
        sentrix_id,
        reference_umap[:-5],
        ENDING["umap_xlsx"],
    )

    url = UMAP_LINK % reference_umap
    request.urlretrieve(url, umap_coordinates_local)

    cnv_local = composite_path(NANODIP_REPORTS, sentrix_id, ENDING["cnv_pdf"])
    url = CNV_LINK % sentrix_id
    request.urlretrieve(url, cnv_local)

    image = convert_from_path(cnv_local)[0]
    image.save(cnv_local.replace("pdf", "png"), "png")


def calculate_std(reference):
    """Calculate sorted standard deviations with GPU (if present) for
    a particular reference dataset and save to disk.
    """
    if not os.path.exists(EPIDIP_TMP):
        os.makedirs(EPIDIP_TMP)

    if gpu_enabled():
        # Get unified pool
        pool = cupy.cuda.MemoryPool(cupy.cuda.memory.malloc_managed)
        # Set unified pool as default allocator
        cupy.cuda.set_allocator(pool.malloc)
        # Release GPU memory
        pool.free_all_blocks()
    else:
        logger.info("Probably no CUDA device present.")

    specimen_bin_files = [
        composite_path(BETA_VALUES, s, ENDING["betas_bin"])
        for s in reference.specimens
    ]
    specimens_cnt = len(reference.specimens)
    cpg_cnt = len(reference.cpg_sites)

    # Determine size of beta value array. Number of CpGs is typically fixed,
    # Number of cases is variable
    block_size = GPU_RAM_USAGE // (GPU_FLOAT_SIZE * len(reference.specimens))

    # TODO BUG: Segmentation fault (core dumped)
    # Reproducing bug: delete all files in tmp/epidip then plot reference MNG
    # afterwards plot GSE.
    # print("**************************0*******************************")
    # print("specimens_cnt=", specimens_cnt, "block_size=", block_size)
    # if reference.name == "GSE90496_IfP01":
        # import pdb; pdb.set_trace()
    # Initialize memory for loop.
    beta_values_xp = xp.full(
        [specimens_cnt, block_size], -1, dtype=float, order="C"
    )
    # print("**************************1*******************************")
    beta_stds_xp = xp.array([])

    # Break data into blocks along the cpg-columns, adjusted to
    # GPU RAM availability.
    for col0 in tqdm(range(0, cpg_cnt, block_size), desc="Calculating stds"):
        col1 = min(cpg_cnt, col0 + block_size)
        # Will be =block_size except for the last run.
        d_col = col1 - col0
        for idx, file_ in enumerate(specimen_bin_files):
            beta_values_xp[idx][:d_col] = xp.fromfile(
                file_, count=d_col, offset=col0, dtype=float
            )
        # Replace nan with 0.49
        beta_values_xp = xp.nan_to_num(beta_values_xp, nan=0.49)
        beta_stds_xp = xp.append(
            beta_stds_xp,
            xp.std(beta_values_xp, axis=0, dtype=float)[:d_col],
        )

    # Convert cupy to numpy array
    if isinstance(beta_stds_xp, cupy.ndarray):
        beta_stds = beta_stds_xp.get()
    else:
        beta_stds = beta_values_xp

    # Standard deviations >1 are useless (typically INF values)
    beta_stds[(beta_stds > 1) | (np.isnan(beta_stds))] = 0
    std_bin = composite_path(EPIDIP_TMP, reference.name, ENDING["stdarr_bin"])
    beta_stds.tofile(std_bin)

    # Create data frame containing cpg sites with stds
    beta_value_df = pd.DataFrame(columns=["index", "cpg_site", "std"])
    beta_value_df.cpg_site = reference.cpg_sites
    beta_value_df["index"] = range(0, cpg_cnt)
    beta_value_df["std"] = beta_stds
    # Sort descending by std
    beta_value_df.sort_values(
        by="std",
        axis=0,
        ascending=False,
        inplace=True,
        kind="quicksort",
        na_position="last",
        ignore_index=False,
        key=None,
    )
    std_sorted = composite_path(
        EPIDIP_TMP, reference.name, ENDING["stdsortarr_bin"]
    )
    beta_value_df.to_csv(path_or_buf=std_sorted, index=False)

    # Need to release GPU memory explicitly
    del beta_stds
    del beta_stds_xp
    del beta_values_xp

    # Release GPU memory
    if gpu_enabled():
        pool.free_all_blocks()

    return beta_value_df


def top_variable_cpgs(reference, nr_top_cpgs):
    """Returns the 'nr_top_cpgs' most variable CpG's of 'reference',
    measured by standard deviation.
    """
    std_sorted = composite_path(
        EPIDIP_TMP, reference.name, ENDING["stdsortarr_bin"]
    )
    # Calculates the stds and stores them in a file if it does not
    # already exist or is older than 1 hour.
    if (
        not os.path.exists(std_sorted)
        or time.time() - os.path.getmtime(std_sorted) > 60 * 60
    ):
        calculate_std(reference)
    beta_value_df = pd.read_csv(std_sorted)
    return beta_value_df["cpg_site"][:nr_top_cpgs]
