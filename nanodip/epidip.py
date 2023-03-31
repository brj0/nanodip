"""
### EpiDiP Functions
"""

# start_external_modules
import logging
import os
import cupy
from tqdm import tqdm
import numpy as np
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
from nanodip.plots import (
    UMAPData,
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


def calculate_std(reference_id):
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
    reference = Reference(reference_id)

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

    # Convert cupy to numpy array
    if isinstance(beta_stds, cupy.ndarray):
        beta_stds = beta_stds.get()

    # Standard deviations >1 are useless (typically INF values)
    beta_stds[beta_stds > 1] = 0
    std_bin = composite_path(EPIDIP_TMP, reference_id, ENDING["stdarr_bin"])
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
    std_sorted_csv = composite_path(
        EPIDIP_TMP, reference_id, ENDING["stdsortarr_bin"]
    )
    beta_value_df.to_csv(path_or_buf=std_sorted_csv, index=False)

    # Need to release GPU memory explicitly
    del beta_value_df
    del beta_stds
    del beta_values

    # Release GPU memory
    if gpu_enabled():
        pool.free_all_blocks()


reference_id = "GSE90496_IfP01"
sample_id = "all_samples"
# calculate_std(reference_id)


umap_data = UMAPData(sample_id, reference_id)
std_sorted_csv = composite_path(
    EPIDIP_TMP, reference_id, ENDING["stdsortarr_bin"]
)
beta_value_df = pd.read_csv(std_sorted_csv)


def epidipUmap(
    referenceName, referenceStdName, topN
):  # calculate UMAP plot from files in a given reference set for topN probes

    if not os.path.exists(epidipTmp):
        os.makedirs(epidipTmp)
    stdFile = epidipTmp + "/" + referenceStdName + "_stdArray.bin"
    stdSortFileCsv = epidipTmp + "/" + referenceStdName + "_stdSortArray.csv"
    epidipUmapXlsx = (
        epidipTmp
        + "/"
        + datetimestringnow()
        + "_EpiDiP_"
        + str(topN)
        + "_"
        + referenceName.replace(".xlsx", "")
        + "_"
        + referenceStdName.replace(".xlsx", "")
        + ".xlsx"
    )
    binSuffix = "_betas_filtered.bin"
    betaStdDf = pandas.read_csv(stdSortFileCsv, header="infer", sep=",")
    referenceString = referenceName.replace(".xlsx", "")
    referenceSheetFile = (
        referenceDir + "/" + referenceName
    )  # load reference annotation
    referenceSheet = pandas.read_excel(
        referenceSheetFile,
        header=None,
        names=["SentrixID", "MethClass", "MethText"],
    )
    binFiles = pandas.DataFrame(
        listdir(binDir)
    )  # collect reference case binary file names
    binFiles.columns = ["binFileName"]  # name first column
    binFiles["SentrixID"] = binFiles.apply(
        lambda row: row.binFileName.replace(binSuffix, ""), axis=1
    )  # get SentrixID with string operation on dataframe
    binFiles["binPath"] = binFiles.apply(
        lambda row: binDir + "/" + row.binFileName, axis=1
    )  # get Path with string operation on dataframe
    referenceSheet = referenceSheet.merge(
        binFiles, on="SentrixID", how="inner"
    )  # get overlap between reference list and available bin files
    numCases = referenceSheet.shape[0]
    floatSize = 8  # size per float in GPU RAM (tested on AGX Xavier)
    betaValues = numpy.full(
        [numCases, topN], -1, dtype="float32", order="C"
    )  # create fixed-size cupy array filled with -1
    float64bytes = 8  # binary float64 representation
    betaStdDf = betaStdDf[0:topN]
    betaStdDf.sort_values(
        by="binIndex",
        axis=0,
        ascending=True,
        inplace=True,
        kind="quicksort",
        na_position="last",
        ignore_index=False,
        key=None,
    )  # sort the topN lines betaStdDf offsets to facilitate faster loading from permanent storage; rewinding is slow
    betaStdDf["binOffset"] = betaStdDf.apply(
        lambda row: row.binIndex * float64bytes, axis=1
    )  # pre-calculate offsets (in bytes)
    ind = list(betaStdDf["binOffset"])
    p_bar0 = tqdm(range(numCases))
    c = 0
    for f in p_bar0:
        with open(referenceSheet["binPath"][f], "rb") as b:
            buf = bytearray()
            for i in ind:
                b.seek(i)  # offset (pre-calculated)
                buf += b.read(float64bytes)  # read bytes into buffer
        betaValues[c] = numpy.float32(numpy.frombuffer(buf, dtype="float64"))
        c += 1
    betaValues = numpy.nan_to_num(
        betaValues, nan=0.49
    )  # replace nan with 0.49
    betaValuesDf = pandas.DataFrame(betaValues)
    embedding = pandas.DataFrame(umap.UMAP().fit_transform(betaValues))
    del betaValues
    del betaStdDf
    embedding.columns = ["UMAP 0", "UMAP 1"]
    referenceSheet = referenceSheet.join(embedding)
    referenceSheet.to_excel(epidipUmapXlsx)
    return epidipUmapXlsx
