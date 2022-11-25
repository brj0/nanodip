"""
### EpiDiP Functions
"""

# start_external_modules
import logging
import os
from tqdm import tqdm
import pandas as pd
# end_external_modules

# start_internal_modules
from nanodip.config import (
    BETA_VALUES,
    ENDING,
    EPIDIP_TMP,
    GPU_FLOAT_SIZE,
    GPU_RAM_USAGE,
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


def gpu_enabled():
    """Tests if CUDA device is present."""
    try:
        import cupy
        cupy.cuda.Device()
        return True
    except:
        return False

# Import cupy if available.
if gpu_enabled():
    import cupy as xp
else:
    import numpy as xp

def calculate_std(reference_name):
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

    # Collect reference case binary file names
    reference = Reference(reference_name)

    specimen_bin_files = [
        composite_path(BETA_VALUES, s, ENDING["betas_bin"])
        for s in reference.specimens
    ]
    specimens_cnt = len(reference.specimens)
    cpg_cnt = len(reference.cpg_sites)

    # Determine size of beta value array. Number of CpGs is typically fixed,
    # Number of cases is variable
    block_size = GPU_RAM_USAGE // (GPU_FLOAT_SIZE * len(reference.specimens))

    # Initialize memory for loop.
    beta_values = xp.full(
        [specimens_cnt, block_size], -1, dtype=float, order="C"
    )
    beta_stds = xp.array([])

    # Break data into blocks along the cpg-columns, adjusted to
    # GPU RAM availability.
    for col0 in tqdm(range(0, cpg_cnt, block_size), desc="Calculating stds"):
        col1 = min(cpg_cnt, col0 + block_size)
        # Will be =block_size except for the last run.
        d_col = col1 - col0
        for idx, file_ in enumerate(specimen_bin_files):
            beta_values[idx][:d_col] = xp.fromfile(
                file_, count=d_col, offset=col0, dtype=float
            )
        # Replace nan with 0.49
        beta_values = xp.nan_to_num(beta_values, nan=0.49)
        beta_stds = xp.append(
            beta_stds,
            xp.std(beta_values, axis=0, dtype=float)[:d_col],
        )

    # Standard deviations >1 are useless (typically INF values)
    beta_stds[beta_stds > 1] = 0
    std_bin = composite_path(EPIDIP_TMP, reference_name, ENDING["stdarr_bin"])
    beta_stds.tofile(std_bin)

    # Create data frame containing cpg sites with stds
    beta_value_df = pd.DataFrame(columns=["index", "cpg_site", "std"])
    beta_value_df.cpg_site = reference.cpg_sites
    beta_value_df.index = range(0, cpg_cnt)
    if isinstance(beta_stds, cupy.ndarray):
        # get is required for cupy arrays only
        beta_value_df.std = beta_stds.get()
    else:
        beta_value_df.std = beta_stds

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
    std_sorted_csv = composite_path(
        EPIDIP_TMP, reference_name, ENDING["stdsortarr_bin"]
    )
    beta_value_df.to_csv(path_or_buf=std_sorted_csv, index=False)

    # Need to release GPU memory explicitly
    del beta_value_df
    del beta_stds
    del beta_values

    # Release GPU memory
    if gpu_enabled():
        pool.free_all_blocks()

reference_name = "GSE90496_IfP01"
calculate_std(reference_name)
