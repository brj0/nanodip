#!/usr/bin/env python

"""
## NanoDiP all-in-one Jupyter Notebook
*J. Hench, S. Frank, and C. Hultschig, Neuropathology, IfP Basel, 2021*

This software is provided free of charge and warranty; by using it you agree
to do this on your own risk. The authors shall not be held liable for any
damage caused by this software. We have assembled this and tested it to the
best of our knowledge.

The purpose of NanoDiP (Nanopore Digital Pathology) is to compare low-coverage
Nanopore sequencing data from natively extracted DNA sequencing runs against
a flexibly adaptable collection of 450K/850K Illumina Infinium Methylation
array data. These data have to be preprocessed into binary beta value files;
this operation is performed in R (uses minfi to read raw array data) and
outputs bindary float files (one per dataset). These beta values files (e.g.,
204949770141_R03C01_betas_filtered.bin) are named according to the array ID
(Sentrix ID) followed by the suffix. A collection of betas_filtered.bin files
can be provided in a static manner and XLSX (Microsoft Excel) tables can be
used to select a subset thereof alongside a user-defined annotation. The
corresponding datasets will be loaded into memory and then serve as the
reference cohort to which the Nanopore data are compared by dimension reduction
(UMAP). This comparison is optimized for speed and low resource consumption so
that it can run on the computer that operates the sequencer. The sequencing run
is initiated through the MinKNOW API by this application. Basecalling and
methylation calling occur as background tasks outside this Jupyter Notebook.
User interaction occurs through a web interface based on CherryPy which has
been tested on Chromium web browser. It is advisable to run it locally, there
are no measures to secure the generated website.

In order to use this application properly please make sure to be somewhat
familiar with Jupyter Notebook. To run the software, press the button called
*restart the kernel, re-run the whole notebook (with dialog)* and confirm
execution. Then, in Chromium Browser, navigate to http://localhost:8080/ and
preferably bookmark this location for convenience. In case of errors, you may
just again click the same button *restart the kernel, re-run the whole notebook
(with dialog)*.
___

### Technical Details

* Tested with Python 3.7.5; 3.8.8 fails to load minknow_api in jupyter
  notebook.
* Verified to run on Ubuntu 18.04/Jetpack on ARMv8 and x86_64 CPUs; not
  tested on Windows and Mac OS. The latter two platforms are unsupported, we
  do not intend to support them.
* **CAUTION**: Requires a *patched* version of minknow api, file
  `[VENV]/lib/python3.7/site-packages/minknow_api/tools/protocols.py`.
  Without the patch, the generated fast5 sequencing data will be unreadable
  with f5c or nanopolish (wrong compression algorithm, which is the default in
  the MinKNOW backend).
"""

# Verify running Python version (should be 3.7.5) and adjust jupyter notebook.
import IPython
import os
from IPython.core.display import display, HTML
# set display witdth to 100%
display(HTML("<style>.container { width:100% !important; }</style>"))
os.system('python --version')

"""
## Multithreading Options
Depending on the number of parallel threads/cores of the underlying hardware,
threading options for multithreaded modules need to be set as
environment-specific parameters. One way to do so is through the *os* module.
"""

# execution-wide multithreading options, set according to your hardware. Jetson
# AGX: suggest "2" needs to be set before importing other modules that query
# these parameters
import os
os.environ["NUMBA_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"

"""
## Modules
This section imports the required modules that should have been installed via
pip. Other package managers have not been tested. To install packages, use the
setup script provided with this software or, alternatively, install them one
by one, ideally in a virtual python environment. Note that the MinKNOW API
requires manual patching after installation with pip.
"""

# python_modules_to_import
# start_external_modules
import argparse
import cherrypy
import datetime
import fnmatch
import logging
from minknow_api.manager import Manager
import minknow_api.statistics_pb2
import minknow_api.device_pb2
from minknow_api.tools import protocols
from numba import jit
import pandas as pd
import psutil
import pysam
import shutil
import socket
import subprocess
import sys
import time
import multiprocessing as mp
# end_external_modules

# start_internal_modules
from config import (
    ANALYSIS_EXCLUSION_PATTERNS,
    ANNOTATIONS,
    BARCODE_NAMES,
    BROWSER_FAVICON,
    CHERRYPY_HOST,
    CHERRYPY_PORT,
    DATA,
    DEBUG_MODE,
    EXCLUDED_FROM_ANALYSIS,
    F5C,
    IMAGES,
    MINIMAP2,
    NANODIP_OUTPUT,
    NANODIP_REPORTS,
    NANODIP_VERSION,
    NEEDED_NUMBER_OF_BASES,
    READS_PER_FILE,
    READ_CPG_RSCRIPT,
    REFERENCE_GENOME_FA,
    REFERENCE_GENOME_MMI,
    RESULT_ENDINGS,
    RSCRIPT,
    SAMTOOLS,
    THIS_HOST,
    VERBOSITY,
)
from data import (
    binary_reference_data_exists,
    ReferenceGenome,
)
from plots import (
    CNVData,
    UMAPData,
)
from utils import (
    convert_html_to_pdf,
    date_time_string_now,
    render_template,
)
# end_internal_modules

"""
# No user editable code below
Do not modify the cells below unless you would like to patch errors or create
something new.

## Sections
1. Generic Functions
2. MinKNOW API Functions
3. CNV Plotter
4. UMAP Methylation Plotter
5. User Interface Functions
6. Report Generator
7. CherryPy Web UI
"""

"""
### 1. Generic Functions
"""

def logpr(v,logstring): # logging funcion that reads verbosity parameter
    if v==1:
        print(str(datetime.datetime.now())+": "+str(logstring))


def get_runs():
    """Return list of run folders from MinKNOW data directory sorted by
    modification time.
    """
    runs = []
    for f in os.listdir(DATA):
        if f not in EXCLUDED_FROM_ANALYSIS:
            file_path = os.path.join(DATA, f)
            mod_time = os.path.getmtime(file_path)
            if os.path.isdir(file_path):
                runs.append([f, mod_time])
    # sort based on modif. date
    runs.sort(key=lambda x: (x[1], x[0]), reverse=True)
    # Remove date after sorting
    return [x[0] for x in runs]

def predominant_barcode(sample_name):
    """Returns the predominante barcode within all fast5 files."""
    fast5_files = []
    for root, _, files in os.walk(os.path.join(DATA, sample_name)):
        fast5_files.extend(
            [os.path.join(root, f) for f in files if f.endswith(".fast5")]
        )
    barcode_hits=[]
    for barcode in BARCODE_NAMES:
        barcode_hits.append(
            len([f for f in fast5_files if barcode in f])
        )
    max_barcode = max(barcode_hits)
    if max_barcode > 1:
        predominant_barcode = BARCODE_NAMES[barcode_hits.index(max_barcode)]
    else:
        predominant_barcode = "undetermined"
    return predominant_barcode

def reference_annotations():
    """Return list of all reference annotation files (MS Excel XLSX format)."""
    annotations = []
    for r in os.listdir(ANNOTATIONS):
        if r.endswith(".xlsx"):
            annotations.append(r)
    return annotations


# TODO del
# write the filename of the UMAP reference for the
def writeReferenceDefinition(sampleId,referenceFile):
    # current run into a text file
    with open(NANODIP_REPORTS+'/'+sampleId+'_selected_reference.txt', 'w') as f:
        f.write(referenceFile)

def write_reference_name(sample_id,reference_name):
    """Write the filename of the UMAP reference for the current run into
    a text file.
    """
    path = os.path.join(
        NANODIP_REPORTS, sample_id + "_selected_reference.txt"
    )
    with open(path, "w") as f:
        f.write(reference_name)


def readReferenceDefinition(sampleId): # read the filename of the UMAP reference for the current sample
    try:
        with open(NANODIP_REPORTS+'/'+sampleId+'_selected_reference.txt', 'r') as f:
            referenceFile=f.read()
    except:
        referenceFile=""
    return referenceFile


def writeRunTmpFile(sampleId,deviceId):
    # current run into a text file
    with open(NANODIP_REPORTS+'/'+sampleId+'_'+deviceId+'_runinfo.tmp', 'a') as f:
        try:
            runId=getActiveRun(deviceId)
        except:
            runId="none"
        ro=getThisRunOutput(deviceId,sampleId,runId)
        readCount=ro[0]
        bascalledBases=ro[1]
        overlapCpGs=getOverlapCpGs(sampleId)
        f.write(str(int(time.time()))+"\t"+
                str(readCount)+"\t"+
                str(bascalledBases)+"\t"+
                str(overlapCpGs)+"\n")


def readRunTmpFile(sampleId):
    print("readRunTmpFile not ready")


def getOverlapCpGs(sampleName):
    methoverlapPath=NANODIP_OUTPUT+"/"+sampleName # collect matching CpGs from sample
    methoverlapTsvFiles=[] # find all *methoverlap.tsv files
    for root, dirnames, filenames in os.walk(methoverlapPath):
        for filename in fnmatch.filter(filenames, '*methoverlap.tsv'):
            methoverlapTsvFiles.append(os.path.join(root, filename))
    methoverlap=[]
    first=True
    for f in methoverlapTsvFiles:
        try: # some fast5 files do not contain any CpGs
            m=pd.read_csv(f, delimiter='\t', header=None, index_col=0)
            if first:
                methoverlap=m
                first=False
            else:
                methoverlap=methoverlap.append(m)
        except:
            logpr(VERBOSITY,"empty file encountered, skipping")
    return len(methoverlap)


def f5cOneFast5(sampleId,analyzeOne=True):
    analyzedCount=0
    thisRunDir=DATA+"/"+sampleId
    pattern = '*.fast5'
    fileList = []
    for dName, sdName, fList in os.walk(thisRunDir): # Walk through directory
        for fileName in fList:
            if fnmatch.fnmatch(fileName, pattern): # Match search string
                fileList.append(os.path.join(dName, fileName))
    calledList=[]
    completedCount=0
    maxBcCount=1 # at least 2 "passed" files (>1) need to be present
    targetBc="undetermined"
    for bc in BARCODE_NAMES:
        thisBc=0
        for f in fileList:
            if bc in f:
                if "_pass_" in f:
                    thisBc+=1
        if thisBc > maxBcCount:
            maxBcCount=thisBc
            targetBc=bc
    f5cAnalysisDir=NANODIP_OUTPUT+"/"+sampleId
    if os.path.exists(f5cAnalysisDir)==False:
        os.mkdir(f5cAnalysisDir)
    thisBcFast5=[]
    thisBcFastq=[]
    for f in fileList:
        if targetBc in f:
            q=f.replace(".fast5","").replace("fast5_pass","fastq_pass")+".fastq"
            if os.path.exists(q): # check if accompanying fastq exists
                thisBcFast5.append(f)
                thisBcFastq.append(q)
                thisBcFileName=f.split("/")
                thisBcFileName=thisBcFileName[len(thisBcFileName)-1].replace(".fast5","") # get name prefix (to be the analysis subdir name later)
                thisAnalysisDir=f5cAnalysisDir+"/"+thisBcFileName
                if os.path.exists(thisAnalysisDir)==False:
                    os.mkdir(thisAnalysisDir)
                target5=thisAnalysisDir+"/"+thisBcFileName+".fast5"
                targetq=thisAnalysisDir+"/"+thisBcFileName+".fastq"
                if os.path.exists(target5)==False:
                    os.symlink(f,target5)             # fast5 symlink
                if os.path.exists(targetq)==False:
                    os.symlink(q,targetq)             #fastq symlink
                if os.path.exists(thisAnalysisDir+"/"+thisBcFileName+"-methoverlapcount.txt")==False:
                    if (analyzeOne==True and analyzedCount==0) or analyzeOne==False:
                        cmd=F5C+" index -t 1 --iop 100 -d "+thisAnalysisDir+" "+targetq
                        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) #index, call methylation and get methylation frequencies
                        p.wait()
                        cmd=MINIMAP2+" -a -x map-ont "+REFERENCE_GENOME_MMI+" "+targetq+" -t 4 | "+SAMTOOLS+" sort -T tmp -o "+thisAnalysisDir+"/"+thisBcFileName+"-reads_sorted.bam"
                        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) # get sorted BAM (4 threads)
                        p.wait()
                        cmd=SAMTOOLS+" index "+thisAnalysisDir+"/"+thisBcFileName+"-reads_sorted.bam"
                        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) # index BAM
                        p.wait()
                        cmd=F5C+" call-methylation -B2000000 -K400 -b "+thisAnalysisDir+"/"+thisBcFileName+"-reads_sorted.bam -g "+REFERENCE_GENOME_FA+" -r "+targetq+" > "+thisAnalysisDir+"/"+thisBcFileName+"-result.tsv"
                        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) # set B to 2 megabases (GPU) and 0.4 kreads
                        p.wait()
                        cmd=F5C+" meth-freq -c 2.5 -s -i "+thisAnalysisDir+"/"+thisBcFileName+"-result.tsv > "+thisAnalysisDir+"/"+thisBcFileName+"-freq.tsv"
                        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                        p.wait()
                        cmd=RSCRIPT+" "+READ_CPG_RSCRIPT+" "+thisAnalysisDir+"/"+thisBcFileName+"-freq.tsv "+ILUMINA_CG_MAP+" "+thisAnalysisDir+"/"+thisBcFileName+"-methoverlap.tsv "+thisAnalysisDir+"/"+thisBcFileName+"-methoverlapcount.txt"
                        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                        p.wait()
                        calledList.append(thisBcFileName)
                        analyzedCount+=1
                else:
                    completedCount+=1
    return "Target = "+targetBc+"<br>Methylation called for "+str(calledList)+". "+str(completedCount+analyzedCount)+"/"+str(len(thisBcFast5))


"""
### 2. MinKNOW API Functions
Check https://github.com/nanoporetech/minknow_api for reference.

The following code requires a patched version of the MinKNOW API, install it
from https://github.com/neuropathbasel/minknow_api.
"""


# Construct a manager using the host + port provided. This is used to connect to
def mkManager():
    return Manager(host=THIS_HOST, port=9501, use_tls=False) # the MinKNOW service trough the MK API.



def listMinionPositions(): # list MinION devices that are currenty connected to the system
    manager = mkManager()
    positions = manager.flow_cell_positions() # Find a list of currently available sequencing positions.
    return(positions)   # User could call {pos.connect()} here to connect to the running MinKNOW instance.


def listMinionExperiments(): # list all current and previous runs in the MinKNOW buffer, lost after MinKNOW restart
    manager=mkManager()
    htmlHost="<b>Host: "+THIS_HOST+"</b><br><table border='1'><tr>"
    positions=manager.flow_cell_positions() # Find a list of currently available sequencing positions.
    htmlPosition=[]
    for p in positions:
        htmlPosinfo="<b>-"+str(p)+"</b><br>"
        connection = p.connect()
        mountedFlowCellID=connection.device.get_flow_cell_info().flow_cell_id # return the flow cell info
        htmlPosinfo=htmlPosinfo+"--mounted flow cell ID: <b>" + mountedFlowCellID +"</b><br>"
        htmlPosinfo=htmlPosinfo+"---"+str(connection.acquisition.current_status())+"<br>" # READY, STARTING, sequencing/mux = PROCESSING, FINISHING; Pause = PROCESSING
        protocols = connection.protocol.list_protocol_runs()
        bufferedRunIds = protocols.run_ids
        for b in bufferedRunIds:
            htmlPosinfo=htmlPosinfo+"--run ID: " + b +"<br>"
            run_info = connection.protocol.get_run_info(run_id=b)
            htmlPosinfo=htmlPosinfo+"---with flow cell ID: " + run_info.flow_cell.flow_cell_id +"<br>"
        htmlPosition.append(htmlPosinfo)
    hierarchy = htmlHost
    for p in htmlPosition:
        hierarchy=hierarchy + "<td valign='top'><tt>"+p+"</tt></td>"
    hierarchy=hierarchy+"</table>"
    return(hierarchy)


def getFlowCellID(thisDeviceId): # determine flow cell ID (if any). Note that some CTCs have an empty ID string.
    mountedFlowCellID="no_flow_cell"
    manager=mkManager()
    positions=manager.flow_cell_positions() # Find a list of currently available sequencing positions.
    for p in positions:
        if thisDeviceId in str(p):
            connection = p.connect()
            mountedFlowCellID=connection.device.get_flow_cell_info().flow_cell_id # return the flow cell info
    return mountedFlowCellID


# This cell starts a run on Mk1b devices and perform several checks concerning
# the run protocol.

# modified from the MinKNOW API on https://github.com/nanoporetech/minknow_api (2021-06)
# created from the sample code at
# https://github.com/nanoporetech/minknow_api/blob/master/python/examples/start_protocol.py
# minknow_api.manager supplies "Manager" a wrapper around MinKNOW's Manager
# gRPC API with utilities for querying sequencing positions + offline
# basecalling tools.
# from minknow_api.manager import Manager

# We need 'find_protocol' to search for the required protocol given a kit +
# product code.
# from minknow_api.tools import protocols
def parse_args():
    """Build and execute a command line argument for starting a protocol.

    Returns:
        Parsed arguments to be used when starting a protocol.
    """
    parser = argparse.ArgumentParser(
        description="""
        Run a sequencing protocol in a running MinKNOW instance.
        """
    )
    parser.add_argument(
        "--host",
        default="localhost",
        help="IP address of the machine running MinKNOW (defaults to localhost)",
    )
    parser.add_argument(
        "--port",
        help="Port to connect to on host (defaults to standard MinKNOW port based on tls setting)",
    )
    parser.add_argument(
        "--no-tls", help="Disable tls connection", default=False, action="store_true"
    )
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")

    parser.add_argument("--sample-id", help="sample ID to set")
    parser.add_argument(
        "--experiment-group",
        "--group-id",
        help="experiment group (aka protocol group ID) to set",
    )
    parser.add_argument(
        "--position",
        help="position on the machine (or MinION serial number) to run the protocol at",
    )
    parser.add_argument(
        "--flow-cell-id",
        metavar="FLOW-CELL-ID",
        help="ID of the flow-cell on which to run the protocol. (specify this or --position)",
    )
    parser.add_argument(
        "--kit",
        required=True,
        help="Sequencing kit used with the flow-cell, eg: SQK-LSK108",
    )
    parser.add_argument(
        "--product-code",
        help="Override the product-code stored on the flow-cell and previously user-specified"
        "product-codes",
    )
    # BASECALL ARGUMENTS
    parser.add_argument(
        "--basecalling",
        action="store_true",
        help="enable base-calling using the default base-calling model",
    )
    parser.add_argument(
        "--basecall-config",
        help="specify the base-calling config and enable base-calling",
    )
    # BARCODING ARGUMENTS
    parser.add_argument(
        "--barcoding", action="store_true", help="protocol uses barcoding",
    )
    parser.add_argument(
        "--barcode-kits",
        nargs="+",
        help="bar-coding expansion kits used in the experiment",
    )
    parser.add_argument(
        "--trim-barcodes", action="store_true", help="enable bar-code trimming",
    )
    parser.add_argument(
        "--barcodes-both-ends",
        action="store_true",
        help="bar-code filtering (both ends of a strand must have a matching barcode)",
    )

    parser.add_argument(
        "--detect-mid-strand-barcodes",
        action="store_true",
        help="bar-code filtering for bar-codes in the middle of a strand",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.0,
        help="read selection based on bar-code accuracy",
    )
    parser.add_argument(
        "--min-score-rear",
        type=float,
        default=0.0,
        help="read selection based on bar-code accuracy",
    )

    parser.add_argument(
        "--min-score-mid",
        type=float,
        default=0.0,
        help="read selection based on bar-code accuracy",
    )
    # ALIGNMENT ARGUMENTS
    parser.add_argument(
        "--alignment-reference",
        help="Specify alignment reference to send to basecaller for live alignment.",
    )
    parser.add_argument(
        "--bed-file", help="Specify bed file to send to basecaller.",
    )
    # Output arguments
    parser.add_argument(
        "--fastq",
        action="store_true",
        help="enables FastQ file output, defaulting to 4000 reads per file",
    )
    parser.add_argument(
        "--fastq-reads-per-file",
        type=int,
        default=4000,
        help="set the number of reads combined into one FastQ file.",
    )
    parser.add_argument(
        "--fast5",
        action="store_true",
        help="enables Fast5 file output, defaulting to 4000 reads per file, this will store raw, "
        "fastq and trace-table data",
    )
    parser.add_argument(
        "--fast5-reads-per-file",
        type=int,
        default=4000,
        help="set the number of reads combined into one Fast5 file.",
    )
    parser.add_argument(
        "--bam",
        action="store_true",
        help="enables BAM file output, defaulting to 4000 reads per file",
    )
    parser.add_argument(
        "--bam-reads-per-file",
        type=int,
        default=4000,
        help="set the number of reads combined into one BAM file.",
    )
    # Read until
    parser.add_argument(
        "--read-until-reference", type=str, help="Reference file to use in read until",
    )
    parser.add_argument(
        "--read-until-bed-file", type=str, help="Bed file to use in read until",
    )
    parser.add_argument(
        "--read-until-filter",
        type=str,
        choices=["deplete", "enrich"],
        help="Filter type to use in read until",
    )
    # Experiment
    parser.add_argument(
        "--experiment-duration",
        type=float,
        default=72,
        help="time spent sequencing (in hours)",
    )
    parser.add_argument(
        "--no-active-channel-selection",
        action="store_true",
        help="allow dynamic selection of channels to select pores for sequencing, "
        "ignored for Flongle flow-cells",
    )
    parser.add_argument(
        "--mux-scan-period",
        type=float,
        default=1.5,
        help="number of hours before a mux scan takes place, enables active-channel-selection, "
        "ignored for Flongle flow-cells",
    )
    parser.add_argument(
        "extra_args",
        metavar="ARGS",
        nargs="*",
        help="Additional arguments passed verbatim to the protocol script",
    )
    args = parser.parse_args()
    # Further argument checks
    # Read until must have a reference and a filter type, if enabled:
    if (
        args.read_until_filter is not None
        or args.read_until_reference is not None
        or args.read_until_bed_file is not None
    ):
        if args.read_until_filter is None:
            print("Unable to specify read until arguments without a filter type.")
            sys.exit(1)

        if args.read_until_reference is None:
            print("Unable to specify read until arguments without a reference type.")
            sys.exit(1)

    if args.bed_file and not args.alignment_reference:
        print("Unable to specify `--bed-file` without `--alignment-reference`.")
        sys.exit(1)

    if (args.barcoding or args.barcode_kits) and not (
        args.basecalling or args.basecall_config
    ):
        print(
            "Unable to specify `--barcoding` or `--barcode-kits` without `--basecalling`."
        )
        sys.exit(1)
    if args.alignment_reference and not (args.basecalling or args.basecall_config):
        print("Unable to specify `--alignment-reference` without `--basecalling`.")

        sys.exit(1)
    if not (args.fast5 or args.fastq):
        print("No output (fast5 or fastq) specified")

    return args

def is_position_selected(position, args):
    """Find if the {position} is selected by command line arguments {args}."""

    # First check for name match:
    if args.position == position.name:
        return True

    # Then verify if the flow cell matches:
    connected_position = position.connect()
    if args.flow_cell_id is not None:
        flow_cell_info = connected_position.device.get_flow_cell_info()
        if (
            flow_cell_info.user_specified_flow_cell_id == args.flow_cell_id
            or flow_cell_info.flow_cell_id == args.flow_cell_id
        ):
            return True

    return False


def startRun():
    """Entrypoint to start protocol example."""
    # Parse arguments to be passed to started protocols:
    run_id=""
    args = parse_args()
    #args = parse_args(minknowApiShellArgumentString.split())

    # Specify --verbose on the command line to get extra details about
    if args.verbose:
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    # Construct a manager using the host + port provided:
    #manager = Manager(host=args.host, port=args.port, use_tls=not args.no_tls)
    manager=mkManager()
    errormessage=""

    # Find which positions we are going to start protocol on:
    positions = manager.flow_cell_positions()
    filtered_positions = list(
        filter(lambda pos: is_position_selected(pos, args), positions)
    )

    # At least one position needs to be selected:
    if not filtered_positions:
        errormessage="No positions selected for protocol - specify `--position` or `--flow-cell-id`"
    else:
        protocol_identifiers = {}
        for pos in filtered_positions:
            # Connect to the sequencing position:
            position_connection = pos.connect()

            # Check if a flowcell is available for sequencing
            flow_cell_info = position_connection.device.get_flow_cell_info()
            if not flow_cell_info.has_flow_cell:
                errormessage="No flow cell present in position "+str(pos)
            else:
                # Select product code:
                if args.product_code:
                    product_code = args.product_code
                else:
                    product_code = flow_cell_info.user_specified_product_code
                    if not product_code:
                        product_code = flow_cell_info.product_code

                # Find the protocol identifier for the required protocol:
                protocol_info = protocols.find_protocol(
                    position_connection,
                    product_code=product_code,
                    kit=args.kit,
                    basecalling=args.basecalling,
                    basecall_config=args.basecall_config,
                    barcoding=args.barcoding,
                    barcoding_kits=args.barcode_kits,
                )

                if not protocol_info:
                    print("Failed to find protocol for position %s" % (pos.name))
                    print("Requested protocol:")
                    print("  product-code: %s" % args.product_code)
                    print("  kit: %s" % args.kit)
                    print("  basecalling: %s" % args.basecalling)
                    print("  basecall_config: %s" % args.basecall_config)
                    print("  barcode-kits: %s" % args.barcode_kits)
                    print("  barcoding: %s" % args.barcoding)
                    errormessage="Protocol build error, consult application log."
                else:
                    # Store the identifier for later:
                    protocol_identifiers[pos.name] = protocol_info.identifier

                    # Start protocol on the requested postitions:
                    print("Starting protocol on %s positions" % len(filtered_positions))
                    for pos in filtered_positions:

                        # Connect to the sequencing position:
                        position_connection = pos.connect()

                        # Find the protocol identifier for the required protocol:
                        protocol_identifier = protocol_identifiers[pos.name]

                        # Now select which arguments to pass to start protocol:
                        print("Starting protocol %s on position %s" % (protocol_identifier, pos.name))

                        # Set up user specified product code if requested:
                        if args.product_code:
                            position_connection.device.set_user_specified_product_code(
                                code=args.product_code
                            )

                        # Build arguments for starting protocol:
                        basecalling_args = None
                        if args.basecalling or args.basecall_config:
                            barcoding_args = None
                            alignment_args = None
                            if args.barcode_kits or args.barcoding:
                                barcoding_args = protocols.BarcodingArgs(
                                    args.barcode_kits,
                                    args.trim_barcodes,
                                    args.barcodes_both_ends,
                                    args.detect_mid_strand_barcodes,
                                    args.min_score,
                                    args.min_score_rear,
                                    args.min_score_mid,
                                )

                            if args.alignment_reference:
                                alignment_args = protocols.AlignmentArgs(
                                    reference_files=[args.alignment_reference], bed_file=args.bed_file,
                                )

                            basecalling_args = protocols.BasecallingArgs(
                                config=args.basecall_config,
                                barcoding=barcoding_args,
                                alignment=alignment_args,
                            )

                        read_until_args = None
                        if args.read_until_filter:
                            read_until_args = protocols.ReadUntilArgs(
                                filter_type=args.read_until_filter,
                                reference_files=[args.read_until_reference],
                                bed_file=args.read_until_bed_file,
                                first_channel=None,  # These default to all channels.
                                last_channel=None,
                            )

                        def build_output_arguments(args, name):
                            if not getattr(args, name):
                                return None
                            return protocols.OutputArgs(
                                reads_per_file=getattr(args, "%s_reads_per_file" % name)
                            )

                        fastq_arguments = build_output_arguments(args, "fastq")
                        fast5_arguments = build_output_arguments(args, "fast5")
                        bam_arguments = build_output_arguments(args, "bam")

                        # print the protocol parameters
                        print("position_connection "+str(position_connection))
                        print("protocol_identifier "+str(protocol_identifier))
                        print("args.sample_id "+str(args.sample_id))
                        print("args.experiment_group "+str(args.experiment_group))
                        print("basecalling_args "+str(basecalling_args))
                        print("read_until_args "+str(read_until_args))
                        print("fastq_arguments "+str(fastq_arguments)) #fastq_arguments OutputArgs(reads_per_file=400)
                        print("fast5_arguments "+str(fast5_arguments)) #fast5_arguments OutputArgs(reads_per_file=400)
                        print("bam_arguments "+str(bam_arguments))
                        print("args.no_active_channel_selection"+str(args.no_active_channel_selection))
                        print("args.mux_scan_period"+str(args.mux_scan_period))
                        print("args.experiment_duration "+str(args.experiment_duration))
                        print("args.extra_args "+str(args.extra_args))  # Any extra args passed.

                        # Now start the protocol:
                        run_id = protocols.start_protocol(
                            position_connection,
                            protocol_identifier,
                            sample_id=args.sample_id,
                            experiment_group=args.experiment_group,
                            basecalling=basecalling_args,
                            read_until=read_until_args,
                            fastq_arguments=fastq_arguments,
                            fast5_arguments=fast5_arguments,
                            bam_arguments=bam_arguments,
                            disable_active_channel_selection=args.no_active_channel_selection,
                            mux_scan_period=args.mux_scan_period,
                            experiment_duration=args.experiment_duration,
                            args=args.extra_args,  # Any extra args passed.
                        )

                        #print("Started protocol %s" % run_id)
    return errormessage+run_id # one of them should be ""


def stopRun(minionId): # stop an existing run (if any) for a MinION device
    manager=mkManager()
    positions = list(manager.flow_cell_positions())
    filtered_positions = list(filter(lambda pos: pos.name == minionId, positions))
    # Connect to the grpc port for the position:
    connection = filtered_positions[0].connect()
    protocols = connection.protocol.list_protocol_runs()
    bufferedRunIds = protocols.run_ids
    thisMessage="No protocol running, nothing was stopped."
    c=0
    for b in bufferedRunIds:
        try:
            connection.protocol.stop_protocol()
            thisMessage="Protocol "+b+" stopped on "+minionId+"."
        except:
            c=c+1
    return thisMessage


# from minknow_api demos, start_seq.py
def is_position_selected(position, args):
    """Find if the {position} is selected by command line arguments {args}."""
    if args.position == position.name: # First check for name match:
        return True
    connected_position = position.connect()  # Then verify if the flow cell matches:
    if args.flow_cell_id is not None:
        flow_cell_info = connected_position.device.get_flow_cell_info()
        if (flow_cell_info.user_specified_flow_cell_id == args.flow_cell_id
            or flow_cell_info.flow_cell_id == args.flow_cell_id):
            return True
    return False


def getMinKnowApiStatus(deviceString): # MinKNOW status per device
    replyString=""
    testHost="localhost"
    manager=mkManager()
    positions = list(manager.flow_cell_positions())
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    # determine if anything is running and the kind of run, via set temperature
    replyString=replyString+"acquisition.get_acquisition_info().state: "+str(connection.acquisition.get_acquisition_info().state)+"<br>"
    replyString=replyString+"acquisition.current_status(): "+str(connection.acquisition.current_status())+"<br>"
    replyString=replyString+"minion_device.get_settings().temperature_target.min: "+str(connection.minion_device.get_settings().temperature_target.min)+"<br>"
    replyString=replyString+"device.get_temperature(): " + str(connection.device.get_temperature().minion.heatsink_temperature)+"<br>"
    replyString=replyString+"device.get_bias_voltage(): " + str(connection.device.get_bias_voltage())+"<br>"
    return replyString



def getActiveRun(deviceString):
    manager=mkManager()
    positions = list(manager.flow_cell_positions())
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    try:
        activeRun=connection.acquisition.get_current_acquisition_run().run_id # error if no acquisition is running, same as with acquisitio.current_status(), no acquisition until temperature reached
    except:
        activeRun="none"
    return activeRun


def getRealDeviceActivity(deviceString):            # seq. runs: 34 degC and flow cell checks 37 degC target
    manager=mkManager()                             # temperatures seem to be the only way to determine if
    positions = list(manager.flow_cell_positions()) # a device has been started
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    targetTemp=str(connection.minion_device.get_settings().temperature_target.min)
    returnValue=""
    if targetTemp=="34.0":
        returnValue="sequencing"
    elif targetTemp=="37.0":
        returnValue="checking flow cell"
    elif targetTemp=="35.0":
        returnValue="idle"
    return returnValue


def getThisRunState(deviceString): # obtain further information about a particular device / run
    manager=mkManager()
    positions = list(manager.flow_cell_positions())
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    try:
        thisRunState="Run state for "+deviceString+": "
        thisRunState=thisRunState+str(connection.protocol.get_current_protocol_run().state)+"/"
        thisRunState=thisRunState+str(connection.acquisition.get_acquisition_info().state)
    except:
        thisRunState="No state information in MinKNOW buffer for "+deviceString
    return thisRunState



def getThisRunSampleID(deviceString): # get SampleID from MinKNOW by device, only available after data
    manager=mkManager()               # acquisition as been initiated by MinKNOW.
    positions = list(manager.flow_cell_positions())
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    try:
        thisRunSampleID=connection.protocol.get_current_protocol_run().user_info.sample_id.value
    except:
        thisRunSampleID="No sampleId information in MinKNOW buffer for "+deviceString
    return thisRunSampleID



def getThisRunYield(deviceString): # get run yield by device. The data of the previous run will remain
    manager=mkManager()            # in the buffer until acquisition (not just a start) of a new run
    positions = list(manager.flow_cell_positions()) # have been initiated.
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    try:
        acqinfo=connection.acquisition.get_acquisition_info()
        thisRunYield="Run yield for "+deviceString+"("+acqinfo.run_id+"):&nbsp;"
        thisRunYield=thisRunYield+str(acqinfo.yield_summary)
    except:
        thisRunYield="No yield information in MinKNOW buffer for "+deviceString
    return thisRunYield



def getThisRunOutput(deviceString,sampleName,runId): # get run yield by device, sampleName, runId
    thisRunOutput=[-1,-1] # defaults in case of error / missing information
    manager=mkManager()            # in the buffer until acquisition (not just a start) of a new run
    positions = list(manager.flow_cell_positions()) # have been initiated.
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    readCount=-3
    calledBases=-3
    if getThisRunSampleID(deviceString)==sampleName: # check that runID and sampleID match
        readCount=-4
        calledBases=-4
        if connection.acquisition.get_current_acquisition_run().run_id==runId:
            if connection.acquisition.current_status()!="status: READY": # i.e., working
                try:
                    acq=connection.acquisition.get_acquisition_info()
                    readCount=acq.yield_summary.basecalled_pass_read_count
                    calledBases=acq.yield_summary.basecalled_pass_bases
                except:
                    readCount=-5
                    calledBases=-5
    thisRunOutput=[readCount,calledBases]
    return thisRunOutput # shall be a list


def getThisRunEstimatedOutput(deviceString,sampleName,runId): # get run yield by device, sampleName, runId
    thisRunOutput=[-1,-1] # defaults in case of error / missing information
    manager=mkManager()            # in the buffer until acquisition (not just a start) of a new run
    positions = list(manager.flow_cell_positions()) # have been initiated.
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    readCount=-3
    calledBases=-3
    if getThisRunSampleID(deviceString)==sampleName: # check that runID and sampleID match
        readCount=-4
        calledBases=-4
        if connection.acquisition.get_current_acquisition_run().run_id==runId:
            if connection.acquisition.current_status()!="status: READY": # i.e., working
                try:
                    acq=connection.acquisition.get_acquisition_info()
                    readCount=acq.yield_summary.basecalled_pass_read_count
                    calledBases=acq.yield_summary.estimated_selected_bases
                except:
                    readCount=-5
                    calledBases=-5
    thisRunOutput=[readCount,calledBases]
    return thisRunOutput # shall be a list


def getThisRunInformation(deviceString): # get current run information. Only available after data acquisition
    manager=mkManager()                  # has started.
    positions = list(manager.flow_cell_positions())
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    try:
        thisRunInfo="Run information for "+deviceString+"<br><br>"+str(connection.protocol.get_current_protocol_run())
    except:
        thisRunInfo="No protocol information in MinKNOW buffer for "+deviceString
    return thisRunInfo


def thisRunWatcherTerminator(deviceString,sampleName):
    realRunId=getActiveRun(deviceString) #
    currentBases=getThisRunEstimatedOutput(deviceString,sampleName,realRunId)[1]
    currentBasesString=str(round(currentBases/1e6,2))
    wantedBasesString=str(round(NEEDED_NUMBER_OF_BASES/1e6,2))
    myString="<html><head>"
    myString=myString+"<title>"+currentBasesString+"/"+wantedBasesString+"MB:"+sampleName+"</title>"
    if currentBases < NEEDED_NUMBER_OF_BASES: # don't refresh after showing the STOP state
        myString=myString+"<meta http-equiv='refresh' content='10'>"
    myString=myString+"</head><body>"
    myString=myString+"<b>Automatic run terminator</b> for sample <b>"+sampleName+ "</b>, run ID="+realRunId+" on "+deviceString+" when reaching "+wantedBasesString+" MB, now "+currentBasesString+" MB"
    myString=myString+"<hr>"
    myString=myString+"Last refresh at "+date_time_string_now()+".<hr>"
    if currentBases > NEEDED_NUMBER_OF_BASES:
        stopRun(deviceString)
        myString=myString+"STOPPED at "+date_time_string_now()
    elif currentBases==0:
        myString=myString+"IDLE / MUX / ETC"
    else:
        myString=myString+"RUNNING"
    myString=myString+"</body></html>"
    return myString


"""
### 3. CNV Plotter
"""




"""
### 4. UMAP Methylation Plotter
"""

"""
### 5. Report Generator
"""


"""
### 6. User Interface Functions
"""


# TODO changed
# String patterns in sample names that exclude data from downstream analysis,
# e.g., test runs
def analysis_launch_table():
    """Presents a html table from which analyses can be started in a post-hoc
    manner.
    """
    analysis_runs = [run for run in get_runs() if
        not any(pattern in run for pattern in ANALYSIS_EXCLUSION_PATTERNS)]
    annotations = reference_annotations()
    table = f"""
        <tt>
        <font size='-2'>
        <table border=1>
        <thead>
        <tr>
            <th align='left'><b>Sample ID </b></th>
            <th align='left'><b>CpGs</b></th>
            <th align='left'><b>CNV</b></th>"""
    for a in annotations:
        table += f"""
            <th align='left'>
                <b>UMAP against<br>{a.replace(".xlsx", "")}</b>
            </th>"""
    table += """
        </tr>
        </thead>
        <tbody>"""
    for _, run in enumerate(analysis_runs):
        table += f"""
        <tr>
            <td>{run}</td>
            <td>
            <a href='./analysisLauncher?functionName=methylationPoller&sampleName={run}&refAnno=None'
            target='_blank' rel='noopener noreferrer' title='{run}: CpGs'>
                get CpGs
            </a>
            </td>
            <td>
            <a href='./analysisLauncher?functionName=cnvplot&sampleName={run}&refAnno=None'
                target='_blank' rel='noopener noreferrer' title='{run}: CNV'>
                    plot CNV
            </a>
            </td>"""
        for a in annotations:
            table += f"""
            <td>
            <a href='./analysisLauncher?functionName=umapplot&sampleName={run}&refAnno={a}'
            target='_blank' rel='noopener noreferrer'
            title='{run}: {a.replace(".xlsx", "")}'>
                plot UMAP
            </a>&nbsp;
            <a href='./makePdf?sampleName={run}&refAnno={a}' target='_blank'
            rel='noopener noreferrer' title='{run}: {a.replace(".xlsx", "")}'>
                make PDF
            </a>
            </td>"""
        table += """
        </tr>"""
    table += """
        </tbody>
        </table>
        </font>
        </tt>"""
    return table

def get_all_results():
    """Return list of all analysis result files in report directory sorted
    by modification time."""
    files = []
    for f in os.listdir(NANODIP_REPORTS):
        for e in RESULT_ENDINGS.values():
            if f.endswith(e):
                mod_time = os.path.getmtime(
                    os.path.join(NANODIP_REPORTS, f)
                )
                files.append([f, mod_time])
    files.sort(key=lambda x: (x[1], x[0]), reverse=True)
    return [f[0] for f in files]

def livePage(deviceString): # generate a live preview of the data analysis with the current PNG figures
    thisSampleID=getThisRunSampleID(deviceString) # if there is a run that produces data, the run ID will exist
    thisSampleRef=readReferenceDefinition(thisSampleID).replace(".xlsx", "")
    cnvPlotPath="reports/"+thisSampleID+"_CNVplot.png"
    umapAllPlotPath="reports/"+thisSampleID+"_"+thisSampleRef+"_UMAP_all.png"
    umapAllPlotlyPath="reports/"+thisSampleID+"_"+thisSampleRef+"_UMAP_all.html"
    umapTopPlotPath="reports/"+thisSampleID+"_"+thisSampleRef+"_UMAP_top.png"
    ht="<html><body><tt>sample ID: "+thisSampleID+" with reference "+thisSampleRef+"</tt><br>"
    ht=ht+"<a href='"+cnvPlotPath+"' target='_blank'><img align='Top' src='"+cnvPlotPath+"' width='50%' alt='CNV plot will appear here'></a>"
    ht=ht+"<a href='"+umapAllPlotlyPath+"' target='_blank'><img align='Top' src='"+umapAllPlotPath+"' width='50%' alt='UMAP plot will appear here'></a>"
    ht=ht+"</tt></table><body></html>"
    return ht

def methcallLivePage(sampleName): # generate a self-refreshing page to invoke methylation calling
    ht="<html><head><title>MethCaller: "+sampleName+"</title>"
    ht=ht+"<meta http-equiv='refresh' content='3'></head><body>"
    ht=ht+"last refresh and console output at "+date_time_string_now()+"<hr>shell output<br><br><tt>"
    #ht=ht+calculateMethylationAndBamFromFast5Fastq(sampleName)
    ht=ht+f5cOneFast5(sampleName,analyzeOne=True)
    ht=ht+"</tt></body></html>"
    return ht

#TODO changed
def menuheader(current_page, autorefresh=0):
    """Generate a universal website header for the UI pages that
    contains a simple main menu.
    """
    menu = {
        "index":[
            "Overview",
            "General system information",
        ],
        "listPositions":[
            "Mk1b Status",
            "Live status of all connected Mk1b devices",
        ],
        "startSequencing":[
            "Start seq.",
            "Start a sequencing run on an idle Mk1b device",
        ],
        "startTestrun":[
            "Start test run",
            "Start a test seq. run on an idle Mk1b device to verify that the previous flow cell wash was successful.",
        ],
        "listExperiments":[
            "Seq. runs",
            "List all buffered runs. Will be purged upon MinKNOW backend restart.",
        ],
        "listRuns":[
            "Results",
            "List all completed analysis results",
        ],
        "analyze":[
            "Analyze",
            "Launch data analyses manually, e.g. for retrospective analysis",
        ],
        "about":[
            "About NanoDiP",
            "Version, etc.",
        ],
    }
    html = f"""
        <html>
        <head>
        <title>
            NanoDiP Version {NANODIP_VERSION}
        </title>"""
    if autorefresh > 0:
        html += f"<meta http-equiv='refresh' content='{autorefresh}'>"
    html += """
        </head>
        <body>
        <table border=0 cellpadding=2>
        <tr>
            <td>
                <img src='img/EpiDiP_Logo_01.png' width='40px' height='40px'>
            </td>"""
    for key, value in menu.items():
        selected_color = "#E0E0E0" if current_page == key else "white"
        html += f"""
            <td bgcolor='{selected_color}'>
                <b>
                <a href='{key}' title='{value[1]}'> {value[0]}
                </a>
                </b>
            </td>"""
    html += f"""
        </tr>
        </table>
        <br>"""
    return html

"""
### 6. CherryPy Web UI
The browser-based user interface is based on CherryPy, which contains an
intergrated web server and serves pages locally. Communication between the
service and browser typically generates static web pages that may or may not
contain automatic self refresh commands. In the case of self-refreshing pages,
the browser will re-request a given page with leads to re-execution of the
respective python functions. The main handles to these function are located in
the Web UI cell below.
"""




class UserInterface(object):
    """The CherryPy Web UI Webserver class defines entrypoints and
    function calls.
    """
    # global variables within the CherryPy Web UI
    cpgQueue = 0 # TODO use mutex instead
    umapQueue = 0
    cnvpQueue = 0
    cnv_lock = mp.Lock()
    umap_lock = mp.Lock()

    @cherrypy.expose
    def index(self):
        total, used, free = shutil.disk_usage(DATA)
        sys_stat = {
            "hostname": socket.gethostname(),
            "disk_total": total // (2**30),
            "disk_used": used // (2**30),
            "disk_free": free // (2**30),
            "memory_free": round(
                psutil.virtual_memory().available * 100
                / psutil.virtual_memory().total
            ),
            "cpu": round(psutil.cpu_percent()),
            "cpgs": UserInterface.cpgQueue,
            "cnvp": len([p for p in mp.active_children() if p.name == "cnv"]),
            "umap": len([p for p in mp.active_children() if p.name == "umap"]),
        }
        return render_template("index.html", sys_stat=sys_stat)

    @cherrypy.expose
    def old(self):
        """Titlepage."""
        html = menuheader('index', 15)
        html += "<tt><b>Computer:</b> "
        html += str(socket.gethostname())
        html += "</tt><br><br>"
        return html

    @cherrypy.expose
    def restart(self):
        cherrypy.engine.restart()
        return render_template("restart.html")

    @cherrypy.expose
    def reset_queue(self, queue_name=""):
        html = menuheader('index', 15)
        if queue_name:
            if queue_name == "cpg":
                UserInterface.cpgQueue = 0
            if queue_name == "umap":
                UserInterface.umapQueue = 0
            if queue_name == "cnvp":
                UserInterface.cnvpQueue = 0
            html += queue_name + " queue reset"
        return html


    @cherrypy.expose
    def listPositions(self):
        myString=menuheader(1,10)
        positions=listMinionPositions()
        for pos in positions:
            n=str(pos.name) # pos.state does not tell much other than that the device is connected with USB ("running")
            myString=myString+"<br><iframe src='DeviceStatusLive?deviceString="+n+"' height='200' width='600' title='"+n+"' border=3></iframe>"
            myString=myString+"<iframe src='AnalysisStatusLive?deviceString="+n+"' height='200' width='600' title='"+n+"' border=3></iframe>"
            myString=myString+"<br><a href='DeviceStatusLive?deviceString="+n+"' target='_blank' title='Click to open device status page in new tab or window'>"+n+"</a>"
            myString=myString+", live state: "+getRealDeviceActivity(n)
            activeRun=getActiveRun(n)
            myString=myString+", active run: "+getActiveRun(n)
            if activeRun!="none":
                myString=myString+" <a href='launchAutoTerminator?sampleName="+getThisRunSampleID(n)+"&deviceString="+n+"' target='_blank'>"
                myString=myString+"<br>Click this link to launch automatic run terminator after"+str(round(NEEDED_NUMBER_OF_BASES/1e6))+" MB.</a>"
                myString=myString+"<br><font color=''#ff0000'><a href='stopSequencing?deviceId="+n+"' title='Clicking this will terminate the current run immediately! Use with care!'>terminate manually</a></font>"
            myString=myString+"<br><br>"
        myString=myString+"</body></html>"
        return myString

    @cherrypy.expose
    def status(self):
        positions = [str(pos.name) for pos in listMinionPositions()]
        print("--------------", positions)
        return render_template(
            "status.html",
            positions=positions,
            mega_bases=NEEDED_NUMBER_OF_BASES // 1e6)

    @cherrypy.expose
    def startSequencing(self,deviceId="",sampleId="",runDuration="",referenceFile=""):
        myString=menuheader(2,0)
        if sampleId:
            if float(runDuration)>=0.1:
                sys.argv = ['',
                            '--host','localhost',
                            '--position',deviceId,
                            '--sample-id',sampleId,
                            '--experiment-group',sampleId,
                            '--experiment-duration',runDuration,
                            '--basecalling',
                            '--fastq',
                            '--fastq-reads-per-file',READS_PER_FILE,
                            '--fast5',
                            '--fast5-reads-per-file',READS_PER_FILE,
                            '--verbose',
                            '--kit','SQK-RBK004',
                            '--barcoding',
                            '--barcode-kits','SQK-RBK004']
                realRunId=startRun()
                writeReferenceDefinition(sampleId,referenceFile)
                myString=myString+"sequencing run started for "+sampleId+" on "+deviceId+" as "+realRunId+" with reference "+referenceFile
                myString=myString+"<hr>"+getThisRunInformation(deviceId)
                myString=myString+"<hr><a href='launchAutoTerminator?sampleName="+sampleId+"&deviceString="+deviceId+"'>"
                myString=myString+"Click this link to launch automatic run terminator after"+str(round(NEEDED_NUMBER_OF_BASES/1e6))+" MB.</a> "
                myString=myString+"If you do not start the run terminator, you will have to terminate the run manually, or it will stop after the predefined time."
        else:
            myString=myString+'''<form action="startSequencing" method="GET">
                Select an idle Mk1b:&nbsp;<select name="deviceId" id="deviceId">'''
            positions=listMinionPositions()
            for pos in positions:
                thisPos=pos.name
                if getRealDeviceActivity(thisPos)=="idle":
                    if getFlowCellID(thisPos)!="":
                        myString=myString+'<option value="'+thisPos+'">'+thisPos+': '+getFlowCellID(thisPos)+'</option>'
            myString=myString+'''
                </select>&nbsp; and enter the sample ID:&nbsp;<input type="text" name="sampleId" />
                &nbsp;for&nbsp;<input type="text" name="runDuration" value="72" />&nbsp;hours.
                &nbsp;Reference set&nbsp;<select name="referenceFile" id="referenceFile">'''
            for ref in reference_annotations():
                myString=myString+'<option value="'+ref+'">'+ref+'</option>'
            myString=myString+'&nbsp;<input type="submit" value="start sequencing now"/></form>'
        return myString



    @cherrypy.expose
    def start(self, device_id="", sample_id="",
              run_duration="", reference_id=""):
        start_now = sample_id and float(run_duration) >= 0.1
        if start_now:
            sys.argv = [
                "",
                "--host", "localhost",
                "--position", device_id,
                "--sample-id", sample_id,
                "--experiment-group", sample_id,
                "--experiment-duration", run_duration,
                "--basecalling",
                "--fastq",
                "--fastq-reads-per-file", READS_PER_FILE,
                "--fast5",
                "--fast5-reads-per-file", READS_PER_FILE,
                "--verbose",
                "--kit", "SQK-RBK004",
                "--barcoding",
                "--barcode-kits", "SQK-RBK004",
            ]
            run_id = startRun()
            write_reference_name(sample_id, reference_id)
            return render_template(
                "start.html",
                start_now=start_now,
                test=False,
                sample_id=sample_id,
                reference_id=reference_id,
                device_id=device_id,
                run_id=run_id,
                mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
                run_info=getThisRunInformation(device_id),
            )
        else:
            positions = [p.name for p in listMinionPositions()]
            idle = [p for p in positions if getRealDeviceActivity(p) == "idle"
                and getFlowCellID(p) != ""]
            return render_template(
                "start.html",
                start_now=start_now,
                test=False,
                idle=idle,
                references=reference_annotations(),
            )


    @cherrypy.expose
    def startTestrun(self,deviceId=""):
        myString=menuheader('startTestrun', 0)
        if deviceId:
            sampleId=date_time_string_now()+"_TestRun_"+getFlowCellID(deviceId)
            sys.argv = ['',
                        '--host','localhost',
                        '--position',deviceId,
                        '--sample-id',sampleId,
                        '--experiment-group',sampleId,
                        '--experiment-duration','0.1',
                        '--basecalling',
                        '--fastq',
                        '--fastq-reads-per-file',READS_PER_FILE,
                        '--fast5',
                        '--fast5-reads-per-file',READS_PER_FILE,
                        '--verbose',
                        '--kit','SQK-RBK004',
                        '--barcoding',
                        '--barcode-kits','SQK-RBK004']
            realRunId=startRun()
            myString=myString+"sequencing run started for "+sampleId+" on "+deviceId+" as "+realRunId
            myString=myString+"<hr>"+getThisRunInformation(deviceId)
        else:
            myString=myString+'''<form action="startTestrun" method="GET">
                Select an idle Mk1b:&nbsp;<select name="deviceId" id="deviceId">'''
            positions=listMinionPositions()
            for pos in positions:
                thisPos=pos.name
                if getRealDeviceActivity(thisPos)=="idle":
                    if getFlowCellID(thisPos)!="":
                        myString=myString+'<option value="'+thisPos+'">'+thisPos+': '+getFlowCellID(thisPos)+'</option>'
            myString=myString+'''
                </select>&nbsp;<input type="submit" value="start test run now (0.1h)"/></form>'''
        return myString

    @cherrypy.expose
    def test_run(self, device_id=""):
        if device_id:
            sample_id = (date_time_string_now() + "_TestRun_"
                + getFlowCellID(device_id))
            sys.argv = [
                "",
                "--host", "localhost",
                "--position", device_id,
                "--sample-id", sample_id,
                "--experiment-group", sample_id,
                "--experiment-duration", "0.1",
                "--basecalling",
                "--fastq",
                "--fastq-reads-per-file", READS_PER_FILE,
                "--fast5",
                "--fast5-reads-per-file", READS_PER_FILE,
                "--verbose",
                "--kit", "SQK-RBK004",
                "--barcoding",
                "--barcode-kits", "SQK-RBK004",
            ]
            run_id = startRun()
            return render_template(
                "start.html",
                start_now=True,
                sample_id=sample_id,
                reference_id="TEST",
                device_id=device_id,
                run_id=run_id,
                mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
                run_info=getThisRunInformation(device_id),
            )
        else:
            positions = [p.name for p in listMinionPositions()]
            idle = [p for p in positions if getRealDeviceActivity(p) == "idle"
                and getFlowCellID(p) != ""]
            return render_template(
                "start.html",
                start_now=False,
                test=True,
                idle=idle,
                references=reference_annotations(),
            )

    @cherrypy.expose
    def stopSequencing(self, deviceId=""):
        myString=menuheader('listPositions', 0)
        myString=myString + stopRun(deviceId)
        myString=myString + "<br><br>Click on any menu item to proceed."
        return myString

    @cherrypy.expose
    def listExperiments(self):
        myString=menuheader('listExperiments', 10)
        myString=myString+"Running and buffered experiments:<br>"
        experiments=listMinionExperiments()
        myString=myString+experiments
        return myString

    @cherrypy.expose
    def list_runs(self):
        status = {}
        mounted_flow_cell_id = {}
        current_status = {}
        flow_cell_id = {}
        buffered_run_ids = {}

        manager = mkManager()
        # Find a list of currently available sequencing positions.
        positions = manager.flow_cell_positions()

        for p in positions:
            connection = p.connect()
            # return the flow cell info
            mounted_flow_cell_id[p] = connection.device.get_flow_cell_info(
                ).flow_cell_id
            # READY, STARTING, sequencing/mux = PROCESSING, FINISHING;
            # Pause = PROCESSING
            current_status[p] = connection.acquisition.current_status()
            protocols = connection.protocol.list_protocol_runs()
            buffered_run_ids[p] = protocols.run_ids
            for b in buffered_run_ids[p]:
                run_info = connection.protocol.get_run_info(run_id=b)
                flow_cell_id[(p, b)] = run_info.flow_cell.flow_cell_id
        return render_template(
            "list_runs.html",
            positions=positions,
            host=CHERRYPY_HOST,
            status=status,
            mounted_flow_cell_id=mounted_flow_cell_id,
            current_status=current_status,
            flow_cell_id=flow_cell_id,
            buffered_run_ids=buffered_run_ids,
        )


    @cherrypy.expose
    def results(self):
        files = get_all_results()
        return render_template("results.html", files=files)

    @cherrypy.expose
    def analyze(self):
        myString=menuheader('analyze',0)
        myString=myString+analysis_launch_table()
        return myString

    @cherrypy.expose
    def analysis(self, func="", samp="", ref="", new="False"):
        if func == "":
            analysis_runs = [run for run in get_runs() if not any(pattern in run
                for pattern in ANALYSIS_EXCLUSION_PATTERNS)]
            annotations = [a.replace(".xlsx", "")
                for a in reference_annotations()]
            return render_template(
                "analysis_start.html",
                analysis_runs=analysis_runs,
                annotations=annotations,
            )
        if func == "cnv":
            genome = ReferenceGenome()
            genes = genome.genes.name.to_list()
            return render_template(
                "analysis_cnv.html",
                sample_name=samp,
                genes=genes,
                new=new,
            )
        if func == "umap":
            return render_template(
                "analysis_umap.html",
                sample_name=samp,
                reference_name=ref,
                new=new,
                first_use = not binary_reference_data_exists(),
            )
        else:
            raise cherrypy.HTTPError(404, "URL not found")

    @cherrypy.expose
    def cnv(self, samp, genes="", new="False"):
        t0=time.time()
        print("NEW**********************",new)
        try:
            cnv_plt_data = CNVData(samp)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")

        def make_plot(cnv_data, lock):
            """Plot function for multiprocessing."""
            lock.acquire()
            if not cnv_data.files_on_disk() or new == "True":
                cnv_data.make_cnv_plot()
            lock.release()

        proc = mp.Process(
            target=make_plot,
            args=(cnv_plt_data, UserInterface.cnv_lock),
            name="cnv",
        )
        proc.start()
        proc.join()
        cnv_plt_data.read_from_disk()
        print("CNV=====================", time.time()-t0)

        return cnv_plt_data.plot_cnv_and_genes([genes])

    @cherrypy.expose
    def umap(self, samp, ref, close_up="", new="False"):
        t0=time.time()
        try:
            umap_data = UMAPData(samp, ref)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")

        def make_plot(plt_data, lock):
            """Plot function for multiprocessing."""
            lock.acquire()
            if not plt_data.files_on_disk() or new == "True":
                plt_data.make_umap_plot()
            lock.release()

        proc = mp.Process(
            target=make_plot,
            args=(umap_data, UserInterface.umap_lock),
            name="umap",
        )
        proc.start()
        proc.join()
        umap_data.read_from_disk()

        if close_up == "True":
            return umap_data.cu_plot_json

        return umap_data.plot_json

    @cherrypy.expose
    def umapplot(self, sampleName=None, refAnno=None):
        html = ""
        if sampleName and refAnno:
            while UserInterface.umapQueue > 0:
                time.sleep(2)
            UserInterface.umapQueue += 1
            reference_name = refAnno.replace(".xlsx", "")
            try:
                make_umap_plot(sampleName, reference_name)
                html_error = ""
            except: #TODO which exception?
                html_error = """
                    <b>
                    <font color='#FF0000'>ERROR OCCURRED, PLEASE RELOAD TAB
                    </font>
                    </b>"""
            html += f"""
                <html>
                <head>
                <title>
                    {sampleName} against {refAnno} at {date_time_string_now()}
                </title>
                <meta http-equiv='refresh' content='1;
                    URL=reports/{sampleName}_{reference_name}_UMAP_all.html'>"
                </head>
                <body>
                {html_error}
                Loading UMAP plot. If it fails,
                <a href='reports/{sampleName}_{reference_name}_UMAP_all.html'>
                    click here to load plot
                </a>.
                </body>
                </html>"""
            UserInterface.umapQueue -= 1
        return html

    @cherrypy.expose # TODO crash if files not on disk
    def make_pdf(self, samp=None, ref=None):
        path = os.path.join(NANODIP_REPORTS, samp + "_cpgcount.txt")
        with open(path, "r") as f:
            overlap_cnt = f.read()

        path = os.path.join(NANODIP_REPORTS, samp + "_alignedreads.txt")
        with open(path, "r") as f:
            read_numbers = f.read()

        cnv_path = os.path.join(NANODIP_REPORTS, samp + "_CNVplot.png") #TODO png
        umap_path = os.path.join(
            NANODIP_REPORTS,
            samp + "_" + ref + "_UMAP_top.png",
        )

        html_report = render_template(
            "pdf_report.html",
            sample_name=samp,
            sys_name=socket.gethostname(),
            date=date_time_string_now(),
            barcode=predominant_barcode(samp),
            reads=read_numbers,
            cpg_overlap_cnt=overlap_cnt,
            reference=ref,
            cnv_path=cnv_path,
            umap_path=umap_path,
        )
        report_name = samp + "_" + ref + "_NanoDiP_report.pdf"
        report_path = os.path.join(NANODIP_REPORTS, report_name)
        convert_html_to_pdf(html_report, report_path)
        raise cherrypy.HTTPRedirect(os.path.join("reports", report_name))

    @cherrypy.expose
    def about(self):
        return render_template("about.html")

    @cherrypy.expose
    def DeviceStatusLive(self,deviceString=""):
        currentFlowCellId=getFlowCellID(deviceString)
        myString="<html><head><title>"+deviceString+": "+currentFlowCellId+"</title>"
        try:
            myString=myString+"<meta http-equiv='refresh' content='2'>"
            if getRealDeviceActivity(deviceString)=="sequencing":
                myString=myString+"<body bgcolor='#00FF00'>"
            else:
                myString=myString+"<body>"
            myString=myString+"<b>"+deviceString+": "+currentFlowCellId+"</b><br><tt>"
            myString=myString+getMinKnowApiStatus(deviceString)
        except:
            myString=myString+"<br>No previous device activity, information will appear as soon as the device has been running once in this session.<br>"
        myString=myString+"Sample ID: "+getThisRunSampleID(deviceString)+"<br>"
        myString=myString+getThisRunState(deviceString)
        myString=myString+"<br>"+getThisRunYield(deviceString)
        myString=myString+"</tt></body></html>"
        return myString

    @cherrypy.expose
    def AnalysisStatusLive(self,deviceString=""):
        myString=""
        if deviceString:
            myString=livePage(deviceString)
        return myString

    @cherrypy.expose
    def analysisLauncher(self,functionName="",sampleName="",refAnno=""):
        if functionName and sampleName and refAnno:
            myString="<html><head><title>"+sampleName+" "+functionName+"</title></head><body>"
            myString=myString+functionName+" launched for "+sampleName+" "
            if refAnno!="None":
                myString=myString+"against "+refAnno
            myString=myString+" at "+date_time_string_now()+". "
            myString=myString+"Frame below will display result upon completion, if this tab/window is kept open."
            if refAnno=="None":
                myString=myString+"<br><iframe src='./"+functionName+"?sampleName="+sampleName+"' height='95%' width='100%' title='"+sampleName+"' border=3></iframe>"
            else:
                myString=myString+"<br><iframe src='./"+functionName+"?sampleName="+sampleName+"&refAnno="+refAnno+"' height='95%' width='100%' title='"+sampleName+"' border=3></iframe>"
        else:
            myString="Nothing to launch. You may close this tab now."
        return myString

    @cherrypy.expose
    def analysisPoller(self,sampleName="",deviceString="",runId=""):
        myString="<html><head>"
        if sampleName and deviceString and runId:
                myString=myString+"<title>Poller: "+sampleName+"/"+deviceString+"/"+runId+"</title>"
                myString=myString+"<meta http-equiv='refresh' content='15'>"
                myString=myString+"<body>"
                myString=myString+"Last refresh for "+sampleName+"/"+deviceString+"/"+runId+" at "+date_time_string_now()
                myString=myString+"</body></html>"
                writeRunTmpFile(sampleName,deviceString)
        return myString

    @cherrypy.expose
    def methylationPoller(self,sampleName=""):
        while UserInterface.cpgQueue>0:
            time.sleep(2)
        UserInterface.cpgQueue+=1
        myString=methcallLivePage(sampleName)
        UserInterface.cpgQueue-=1
        return myString

    @cherrypy.expose
    def launchAutoTerminator(self,sampleName="",deviceString=""):
        myString="ERROR"
        if sampleName and deviceString:
            myString=thisRunWatcherTerminator(deviceString,sampleName)
        return myString

def main():
    # Start CherryPy Webserver
    if DEBUG_MODE:
        #set access logging
        cherrypy.log.screen = True
        cherrypy.config.update({'log.screen': True})
    else:
        #set access logging
        cherrypy.log.screen = False
        cherrypy.config.update({'log.screen': False})
        cherrypy.config.update({ "environment": "embedded" })

    print(f"NanoDiP server running at http://{CHERRYPY_HOST}:{CHERRYPY_PORT}")

    cherrypy_config = {
        '/favicon.ico': {
            'tools.staticfile.on': True,
            'tools.staticfile.filename': BROWSER_FAVICON,
        },
        '/img': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': IMAGES,
        },
        '/reports': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': NANODIP_REPORTS,
        },
        '/static': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': os.path.join(os.getcwd(), "static"),
        },
    }
    cherrypy.quickstart(UserInterface(), "/", cherrypy_config)

if __name__ == "__main__":
    main()

"""
### ^^^ LIVE LOG ABOVE ^^^
All CherryPy access will be logged here, including live progress bars for
computationally intense analyses. Detailed access logging is turned off by
default (accessLogging is False), but can be turned on,e.g., for debugging,
in the configuration section at the beginning of this notebook. While it is not
required to have at look at these during normal operation, information
contained in the log may be helpful in troubleshooting. Line numbers in error
messages indicated here typically match those given in the respective Jupyter
Notebook cells.

To preseve these messages, halt the Python kernel, save and close the notebook
to send it for support. This makes sure that the code as well as the error
messages will be preserved.

To launch the user interface, wait until you see a pink log entry that the web
server has started, then navigate to http://localhost:8080.
"""
