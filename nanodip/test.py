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
    BARCODE_NAMES,
    BROWSER_FAVICON,
    CHERRYPY_HOST,
    CHERRYPY_PORT,
    CNV_GRID,
    CNV_LINK,
    CNV_URL_PREFIX,
    CNV_URL_SUFFIX,
    DATA,
    DEBUG_MODE,
    EPIDIP_SERVER,
    EPIDIP_UMAP_COORDINATE_FILES,
    ENDINGS,
    EXCLUDED_FROM_ANALYSIS,
    F5C,
    ILUMINA_CG_MAP,
    IMAGES,
    MINIMAP2,
    NANODIP_OUTPUT,
    NANODIP_REPORTS,
    NANODIP_VERSION,
    NEEDED_NUMBER_OF_BASES,
    PLOTLY_RENDER_MODE,
    READS_PER_FILE,
    REFERENCE_GENOME_FA,
    REFERENCE_GENOME_MMI,
    RESULT_ENDINGS,
    SAMTOOLS,
    THIS_HOST,
    UMAP_PLOT_TOP_MATCHES,
    VERBOSITY,
)
from utils import (
    date_time_string_now,
    discrete_colors,
    composite_path,
)
from data import (
    get_sample_methylation,
    SampleData,
    ReferenceData,
    ReferenceGenome,
    get_reference_methylation,
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
    device_status,
    active_run,
    called_bases,
    run_sample_id,
    start_run,
)
from api import (
    connection_from_device_id,
    parse_args,
    is_position_selected,
    methylation_caller,
)
import config, data, plots, nanodip, utils, api

sample_name = "B2021_48459_20211112_BC10"
sample_name = "B2021_48700_20211112_BC11"
sample_name = "test28"
reference_name = "GSE90496_IfP01"
reference_name = "AllIDATv2_20210804"

sample = SampleData(sample_name)
reference = ReferenceData(reference_name)
# genome = ReferenceGenome()

# data.make_binary_reference_data()
# cnv = CNVData(sample_name)

umapp = UMAPData(sample_name, reference_name)
umapp.read_from_disk()

umapp.sample = SampleData(sample_name)
umapp.reference = ReferenceData(reference_name)
umapp.sample.set_cpg_overlap(reference)
# plot = plots.umap_plot_from_data(
#    sample,
#    reference,
#    umapp.umap_df,
#    close_up=False,
# )
## Write UMAP plot to disk.
# file_path = os.path.join(
#    NANODIP_REPORTS, "__UMAP_TEST__.html"
# )
# plot.write_html(file_path, config=dict({"scrollZoom": True}))


# pie_chart(umapp)


# device_id = "MN26636"
# sample_id = "test"
# run_duration = "0.1"
# start_voltage = "-180"
# run_id = active_run(device_id)
# run_ids=start_run(device_id = "MN26636", sample_id = "test", run_duration = "0.1", start_voltage = "-180")

# l = Devices()
# l.get("a")
# l.get("b")
# l.get("c")
# l.get("f")
# print("l=",l)
# l.pop("c")
# print("l=",l)

from urllib import request

sentrix_id = "201869680197_R07C01"
# wget -O umapLocalFile "http://s1665.rootserver.io/umap_links/UMAP_all_bVals_top_25000.xlsx"
umap_coordinates = "UMAP_all_bVals_top_25000.xlsx"


umap_coordinates_local = composite_path(
    NANODIP_REPORTS, sentrix_id, umap_coordinates[:-5], ENDINGS["umap_xlsx"]
)

URL = EPIDIP_SERVER + umap_coordinates
response = request.urlretrieve(URL, umap_coordinates_local)

cnv_local = composite_path(NANODIP_REPORTS, sentrix_id, ENDINGS["cnv_pdf"])
URL = CNV_LINK % sentrix_id
response = request.urlretrieve(URL, cnv_local)


def downloadEpidipCoordinates(umapFile,sentrixid,referenceFile):
    umapLocalFile=nanodipReportDir+"/"+sentrixid+"_"+referenceFile.replace(".xlsx","")+"_UMAP.xlsx"
    cnvLocalFile=nanodipReportDir+"/"+sentrixid+"_CNVplot.pdf"
    p=subprocess.run(["wget", "-O",
                      umapLocalFile,
                      epidipUmapCoordUrlRoot+umapFile],
                      capture_output=True)
    epidipDownloadStatus='exit status: '+str(p.returncode)+'\nstdout: '+str(p.stdout.decode())+'\nstderr: '+str(p.stderr.decode())
    p=subprocess.run(["wget", "-O",
                      cnvLocalFile,
                      cnvLinkPrefix+sentrixid+cnvLinkSuffix],
                      capture_output=True)
    epidipDownloadStatus=epidipDownloadStatus+'<hr>exit status: '+str(p.returncode)+'\nstdout: '+str(p.stdout.decode())+'\nstderr: '+str(p.stderr.decode())
    p=subprocess.run(["pdftoppm","-png","-f","1","-l","1",
                      cnvLocalFile,
                      cnvLocalFile.replace(".pdf","")],
                      capture_output=True)
    epidipDownloadStatus=epidipDownloadStatus+'<hr>exit status: '+str(p.returncode)+'\nstdout: '+str(p.stdout.decode())+'\nstderr: '+str(p.stderr.decode())
    p=subprocess.run(["mv",
                      cnvLocalFile.replace(".pdf","-1.png"),
                      cnvLocalFile.replace(".pdf",".png")],
                      capture_output=True)
    epidipDownloadStatus=epidipDownloadStatus+'<hr>exit status: '+str(p.returncode)+'\nstdout: '+str(p.stdout.decode())+'\nstderr: '+str(p.stderr.decode())
    return str(epidipDownloadStatus)


