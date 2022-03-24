#!/usr/bin/env python

"""
## NanoDiP all-in-one Jupyter Notebook
*J. Hench, S. Frank, C. Hultschig and J. Brugger, Neuropathology, IfP Basel,
2021*

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
import json
import logging
import grpc
from minknow_api.manager import Manager
import minknow_api.statistics_pb2
import minknow_api.device_pb2
from minknow_api.tools import protocols
from minknow_api.acquisition_pb2 import READY
from numba import jit
import pandas as pd
import psutil
import pysam
import shutil
import socket
import subprocess
import sys
import time
import threading
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
    files_by_ending,
    ReferenceGenome,
)
from plots import (
    CNVData,
    UMAPData,
)
from utils import (
    convert_html_to_pdf,
    date_time_string_now,
    extract_referenced_cpgs,
    render_template,
    url_for,
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
    fast5_files = files_by_ending(DATA, sample_name, ending=".fast5")
    # TODO @HEJU im Original fehlt diese Zeile:
    pass_fast5_files = [f for f in fast5_files if "_pass_" in f]
    barcode_hits=[]
    for barcode in BARCODE_NAMES:
        barcode_hits.append(
            len([f for f in pass_fast5_files if barcode in f])
        )
    max_barcode_cnt = max(barcode_hits)
    if max_barcode_cnt > 1:
        predominant_barcode = BARCODE_NAMES[
            barcode_hits.index(max_barcode_cnt)
        ]
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


def write_reference_name(sample_id, reference_name):
    """Write the filename of the UMAP reference for the current run into
    a text file.
    """
    path = os.path.join(
        NANODIP_REPORTS, sample_id + "_selected_reference.txt"
    )
    with open(path, "w") as f:
        f.write(reference_name)

def read_reference(sample_id):
    """Read the filename of the UMAP reference for the current sample."""
    path = os.path.join(
        NANODIP_REPORTS, sample_id + "_selected_reference.txt"
    )
    try:
        with open(path, "r") as f:
            reference = f.read()
    except FileNotFoundError:
        reference = ""
    return reference

def single_file_methylation_caller(analysis_dir, file_name):
    """Invokes f5c methylation caller on a single fast5/fastq file and
    calculates methylation frequencies and CpG overlaps.

    Args:
        analysis_dir: directory containing fast5 and fastq files.
        file_name: run id.
    """
    base_path = os.path.join(analysis_dir, file_name)
    # Create index file for f5c.
    f5c_index = [
        F5C, "index",
        "-t", "1",
        "--iop", "100",
        "--directory", analysis_dir,
        base_path + ".fastq",
    ]
    # Aligns reads to reference genome and sorts resulting bam files
    # (4 threads).
    seq_align = [
        MINIMAP2,
        "-a",
        "-x", "map-ont",
        REFERENCE_GENOME_MMI,
        base_path + ".fastq",
        "-t", "4",
        "|",
        SAMTOOLS, "sort",
        "-T", "tmp",
        "-o", base_path + "-reads_sorted.bam",
    ]
    # Make bam index for samtools.
    bam_index = [
        SAMTOOLS, "index",
        base_path + "-reads_sorted.bam",
    ]
    # Methylation caller.
    methyl_calling = [
        F5C, "call-methylation",
        #"--disable-cuda=yes",   # TODO for debugging on CPU only. Must be del.
        "-B2000000", "-K400",   # set B to 2 megabases (GPU) and 0.4 kreads
        "-b", base_path + "-reads_sorted.bam",
        "-g", REFERENCE_GENOME_FA,
        "-r", base_path + ".fastq",
        ">", base_path + "-result.tsv",
    ]
    # Calculate methylation frequencies.
    methyl_frequency = [
        F5C, "meth-freq",
        "-c", "2.5",
        "-s",
        "-i", base_path + "-result.tsv",
        ">",
        base_path + "-freq.tsv",
    ]
    commands = [
        f5c_index,
        seq_align,
        bam_index,
        methyl_calling,
        methyl_frequency,
    ]
    # TODO check for successful termination
    for cmd in commands:
        cmd_str = " ".join(cmd)
        p = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE)
        p.wait()
    # Calculate CpG overlap.
    extract_referenced_cpgs(
        base_path + "-freq.tsv",
        base_path + "-methoverlap.tsv",
        base_path + "-methoverlapcount.txt",
    )
    # Write empty textfile to signal successful completion.
    with open(os.path.join(analysis_dir, "done.txt"), "w"):
        pass
    print(f"Methylation calling on {file_name} done.")

def methylation_caller(sample_name, analyze_one=True):
    """Searches for callable fast5/fastq files that have not yet been
    called and invokes methylation calling. Results will be added to
    the NANODIP_OUTPUT directory.

    Args:
        sample_name: Name of sample to be analyzed.
        analyse_one: If True only first fast5/fastq file found
                     will be analyzed.
    """
    # At least 2 "passed" files need to be present.
    barcode = predominant_barcode(sample_name)
    fast5_files = [
        f for f in files_by_ending(DATA, sample_name, ending=".fast5")
        if barcode in f
    ]
    # Analyse in alphanumeric ordering for improved debugging.
    fast5_files.sort()
    def from_5_to_q(fn):
        return fn.replace(
            ".fast5", ".fastq"
        ).replace("fast5_pass", "fastq_pass")

    # TODO @HEJU Daten mit _fail_ werden so nicht ausgeschlossen. Gewollt?
    fast5q_file_pairs = [
        [f, from_5_to_q(f)] for f in fast5_files
        if os.path.exists(from_5_to_q(f))
    ]

    f5c_analysis_dir = os.path.join(NANODIP_OUTPUT, sample_name)
    if not os.path.exists(f5c_analysis_dir):
        os.mkdir(f5c_analysis_dir)

    prev_called = []
    curr_called = []
    not_called = []

    for f5, fq in fast5q_file_pairs:
        file_name = os.path.basename(f5).split(".")[0]
        analysis_dir = os.path.join(f5c_analysis_dir, file_name)
        symlink5 = os.path.join(analysis_dir, file_name + ".fast5")
        symlinkq = os.path.join(analysis_dir, file_name + ".fastq")
        if not os.path.exists(analysis_dir):
            os.mkdir(analysis_dir)
        if not os.path.exists(symlink5):
            os.symlink(f5, symlink5)
        if not os.path.exists(symlinkq):
            os.symlink(fq, symlinkq)
        done = os.path.join(analysis_dir, "done.txt")
        if os.path.exists(done):
            prev_called.append(file_name)
        else:
            not_called.append(
                [analysis_dir, file_name]
            )
    for directory, file_name in not_called:
        single_file_methylation_caller(directory, file_name)
        curr_called.append(file_name)
        if analyze_one:
            break
    num_completed = len(prev_called) + len(curr_called)
    num_fastq = len(fast5q_file_pairs)
    no_callable_left = num_fastq == num_completed
    return {
        "barcode": barcode,
        "called": curr_called,
        "num_completed": num_completed,
        "num_fastq": num_fastq,
        "no_callable_left": no_callable_left,
        "time": date_time_string_now(),
    }

"""
### MinKNOW API Functions
Check https://github.com/nanoporetech/minknow_api for reference.

The following code requires a patched version of the MinKNOW API, install it
from https://github.com/neuropathbasel/minknow_api.
"""

def mk_manager():
    """Construct a manager using the host and port provided. This is
    used to connect to the MinKNOW service trough the MK API.

    minknow_api.manager.Manager:  a wrapper around MinKNOW's Manager
        gRPC API with utilities for querying sequencing positions and
        offline basecalling tools.
    """
    return Manager(host=THIS_HOST, port=9501, use_tls=False)

def minion_positions():
    """Return MinION devices that are currenty connected to the system."""
    manager = mk_manager()
    # Find a list of currently available sequencing positions.
    positions = manager.flow_cell_positions()
    # User could call {posisions.connect()} here to connect to the
    # running MinKNOW instance.
    return positions

def connection_from_device_id(device_id):
    """Returns minion position """
    position = next(
        (pos for pos in minion_positions() if pos.name == device_id),
        False,
    )
    if not position:
        raise ValueError(f"'{device_id}' is not a valid Minion position.")
    connection = position.connect()
    return connection

def flow_cell_id(device_id):
    """Return flow cell ID (if any). Note that some CTCs have an
    empty ID string.
    """
    flow_cell_id = "no_flow_cell"
    positions = minion_positions()
    for p in positions:
        if device_id in p.name:
            connection = p.connect()
            flow_cell_id = connection.device.get_flow_cell_info().flow_cell_id
    return flow_cell_id

def parse_args(args_list):
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
        "--no-tls",
        help="Disable tls connection",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    parser.add_argument(
        "--sample-id",
        help="sample ID to set",
    )
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
    # Basecalling arguments
    parser.add_argument(
        "--basecalling",
        action="store_true",
        help="enable base-calling using the default base-calling model",
    )
    parser.add_argument(
        "--basecall-config",
        help="specify the base-calling config and enable base-calling",
    )
    # Barcoding arguments
    parser.add_argument(
        "--barcoding",
        action="store_true",
        help="protocol uses barcoding",
    )
    parser.add_argument(
        "--barcode-kits",
        nargs="+",
        help="bar-coding expansion kits used in the experiment",
    )
    parser.add_argument(
        "--trim-barcodes",
        action="store_true",
        help="enable bar-code trimming",
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
    # Alignment arguments
    parser.add_argument(
        "--alignment-reference",
        help="Specify alignment reference to send to basecaller for live alignment.",
    )
    parser.add_argument(
        "--bed-file",
        help="Specify bed file to send to basecaller.",
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
    # Read until arguments
    parser.add_argument(
        "--read-until-reference",
        type=str,
        help="Reference file to use in read until",
    )
    parser.add_argument(
        "--read-until-bed-file",
        type=str,
        help="Bed file to use in read until",
    )
    parser.add_argument(
        "--read-until-filter",
        type=str,
        choices=["deplete", "enrich"],
        help="Filter type to use in read until",
    )
    # Experiment arguments
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
    args = parser.parse_args(args_list)
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
        print("Unable to specify '--bed-file' without '--alignment-reference'.")
        sys.exit(1)

    if (args.barcoding or args.barcode_kits) and not (
        args.basecalling or args.basecall_config
    ):
        print(
            "Unable to specify '--barcoding' or '--barcode-kits' without '--basecalling'."
        )
        sys.exit(1)
    if args.alignment_reference and not (args.basecalling or args.basecall_config):
        print("Unable to specify '--alignment-reference' without '--basecalling'.")

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

def start_run(
    device_id="",
    sample_id="",
    run_duration="",
    start_voltage="",
):
    """Start a run on Mk1b devices and perform several checks concerning
    the run protocol.

    Code modified from the MinKNOW API on
    https://github.com/nanoporetech/minknow_api
    (2022-03) created from the sample code at
    https://github.com/nanoporetech/minknow_api/blob/master/python/minknow_api/examples/start_protocol.py

    We need 'find_protocol' to search for the required protocol given a kit
    and product code.
    """
    args_list = [
        "--host",
        "localhost",
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
        "--",  # Required for so-called extra-arguments.
        "--start_bias_voltage", start_voltage,
    ]

    # Parse arguments to be passed to started protocols:
    args = parse_args(args_list)

    # Specify --verbose on the command line to get extra details.
    if args.verbose:
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    # Find which positions we are going to start protocol on:
    positions = [
        pos for pos in minion_positions() if is_position_selected(pos, args)
    ]

    # At least one position needs to be selected:
    if not positions:
        print(
            "No positions selected for protocol - specify "
            "'--position' or '--flow-cell-id'"
        )
        return []

    # Start protocol on the requested postitions:
    print("Starting protocol on %s positions." % len(positions))
    run_ids = []

    for pos in positions:
        # Connect to the sequencing position:
        connection = pos.connect()

        # Check if a flowcell is available for sequencing
        flow_cell_info = connection.device.get_flow_cell_info()
        if not flow_cell_info.has_flow_cell:
            print(f"No flow cell present in position {pos}")
            return []

        # Select product code:
        if args.product_code:
            product_code = args.product_code
        else:
            product_code = flow_cell_info.user_specified_product_code
            if not product_code:
                product_code = flow_cell_info.product_code

        # Find the protocol identifier for the required protocol:
        protocol_info = protocols.find_protocol(
            connection,
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
            print("Protocol build error, consult application log.")
            return []

        # Store the identifier for later:
        protocol_id = protocol_info.identifier

        # Now select which arguments to pass to start protocol:
        print("Starting protocol %s on position %s" % (protocol_id, pos.name))

        # Set up user specified product code if requested:
        if args.product_code:
            connection.device.set_user_specified_product_code(
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
                    reference_files=[args.alignment_reference],
                    bed_file=args.bed_file,
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
        print("connection {connection}")
        print("protocol_id {protocol_id}")
        print("args.sample_id {args.sample_id}")
        print("args.experiment_group {args.experiment_group}")
        print("basecalling_args {basecalling_args}")
        print("read_until_args {read_until_args}")
        print(
            "fastq_arguments {fastq_arguments}"
        )  # fastq_arguments OutputArgs(reads_per_file=400)
        print(
            "fast5_arguments {fast5_arguments}"
        )  # fast5_arguments OutputArgs(reads_per_file=400)
        print("bam_arguments {bam_arguments}")
        print(
            "args.no_active_channel_selection "
            "{args.no_active_channel_selection}"
        )
        print("args.mux_scan_period {args.mux_scan_period}")
        print("args.experiment_duration {args.experiment_duration}")
        print(
            "args.extra_args {args.extra_args}"
        )  # Any extra args passed.

        # Now start the protocol:
        run_id = protocols.start_protocol(
            connection,
            protocol_id,
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
        run_ids.append(run_id)
    return run_ids

def stop_run(device_id):
    """Stop an existing run (if any) for a MinION device and return the
    protocol ID.
    """
    connection = connection_from_device_id(device_id)
    protocols = connection.protocol.list_protocol_runs()
    protocol_id = protocols.run_ids[-1]
    # TODO @HEJU in stopRun wird über bufferedRunIds geloopt. Notwendig?
    try:
        connection.protocol.stop_protocol()
        return protocol_id
    except grpc._channel._InactiveRpcError:
        return None

def is_position_selected(position, args):
    """Find if the {position} is selected by command line arguments
    {args}.
    Function from minknow_api demos, start_seq.py
    """
    if args.position == position.name: # First check for name match:
        return True
    connected_position = position.connect()  # Then verify if the flow cell matches:
    if args.flow_cell_id is not None:
        flow_cell_info = connected_position.device.get_flow_cell_info()
        if (flow_cell_info.user_specified_flow_cell_id == args.flow_cell_id
            or flow_cell_info.flow_cell_id == args.flow_cell_id):
            return True
    return False


def device_status(device_id):
    """MinKNOW status for device."""
    connection = connection_from_device_id(device_id)
    status = {
        "acquisition.get_acquisition_info().state": str(
            connection.acquisition.get_acquisition_info().state
        ),
        "acquisition.current_status()": str(
            connection.acquisition.current_status()
        ),
        "minion_device.get_settings().temperature_target.min": str(
            connection.minion_device.get_settings().temperature_target.min
        ),
        "device.get_temperature()": str(
            connection.device.get_temperature().minion.heatsink_temperature
        ),
        "device.get_bias_voltage()": str(connection.device.get_bias_voltage()),
    }
    return status

def active_run(device_id):
    """Returns active run id."""
    connection = connection_from_device_id(device_id)
    try:
        # Error if no acquisition is running, same as with
        # acquisitio.current_status(), no acquisition until
        # temperature reached
        active_run = connection.acquisition.get_current_acquisition_run().run_id
    except grpc._channel._InactiveRpcError:
        active_run = "none"
    return active_run

def real_device_activity(device_id):
    """Returns device activity by checking the target temperature."""
    connection = connection_from_device_id(device_id)
    target_temp =str(
        connection.minion_device.get_settings().temperature_target.min
    )
    device_activity = {
        "34.0": "sequencing",
        "35.0": "idle",
        "37.0": "checking flow cell",
    }
    return device_activity.get(target_temp, "")

def run_state(device_id):
    """Obtain further information about a particular device / run."""
    connection = connection_from_device_id(device_id)
    try:
        state = f"Run state for {device_id}: "
        state += str(connection.protocol.get_current_protocol_run().state)
        state += "/"
        state += str(connection.acquisition.get_acquisition_info().state)
    except grpc._channel._InactiveRpcError:
        state = f"No state information in MinKNOW buffer for {device_id}"
    return state

def run_sample_id(device_id):
    """Get SampleID from MinKNOW by device, only available after data
    acquisition as been initiated by MinKNOW.
    """
    connection = connection_from_device_id(device_id)
    try:
        sample_id = (
            connection.protocol.get_current_protocol_run().user_info.sample_id.value
        )
    except grpc._channel._InactiveRpcError:
        sample_id = (
            f"No sampleId information in MinKNOW buffer for {device_id}"
        )
    return sample_id

def run_yield(device_id):
    """Get run yield by device. The data of the previous run will remain in
    the buffer until acquisition (not just a start) of a new run have been
    initiated.
    """
    connection = connection_from_device_id(device_id)
    try:
        acq_info = connection.acquisition.get_acquisition_info()
        yield_ = f"Run yield for {device_id} ({acq_info.run_id}):&nbsp;"
        yield_ += str(acq_info.yield_summary)
    # TODO this exception has not been tested.
    except grpc._channel._InactiveRpcError:
        yield_ = f"No yield information in MinKNOW buffer for {device_id}"
    return yield_

def getThisRunEstimatedOutput(deviceString,sampleName,runId): # get run yield by device, sampleName, runId
    thisRunOutput=[-1,-1] # defaults in case of error / missing information
    manager=mk_manager()            # in the buffer until acquisition (not just a start) of a new run
    positions = list(manager.flow_cell_positions()) # have been initiated.
    filtered_positions = list(filter(lambda pos: pos.name == deviceString, positions))
    connection = filtered_positions[0].connect() # Connect to the grpc port for the position
    readCount=-3
    calledBases=-3
    if run_sample_id(deviceString)==sampleName: # check that runID and sampleID match
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

def called_bases(device_id):
    """Returns number of called bases."""
    connection = connection_from_device_id(device_id)
    # Check if device is working.
    if connection.acquisition.current_status().status == READY:
        return 0
    else:
        acquisition = connection.acquisition.get_acquisition_info()
        called_bases = acquisition.yield_summary.estimated_selected_bases
        return called_bases


def run_information(device_id):
    """Get current run information. Only available after data acquisition
    has started.
    """
    connection = connection_from_device_id(device_id)
    try:
        info = (f"Run information for {device_id}<br><br>"
               + str(connection.protocol.get_current_protocol_run()))
    except grpc.RpcError:
        info = f"No protocol information in MinKNOW buffer for {device_id}"
    return info

def thisRunWatcherTerminator(deviceString,sampleName):
    realRunId=active_run(deviceString) #
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
        stop_run(deviceString)
        myString=myString+"STOPPED at "+date_time_string_now()
    elif currentBases==0:
        myString=myString+"IDLE / MUX / ETC"
    else:
        myString=myString+"RUNNING"
    myString=myString+"</body></html>"
    return myString

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
    thisSampleID=run_sample_id(deviceString) # if there is a run that produces data, the run ID will exist
    thisSampleRef=read_reference(thisSampleID).replace(".xlsx", "")
    cnvPlotPath="reports/"+thisSampleID+"_CNVplot.png"
    umapAllPlotPath="reports/"+thisSampleID+"_"+thisSampleRef+"_UMAP_all.png"
    umapAllPlotlyPath="reports/"+thisSampleID+"_"+thisSampleRef+"_UMAP_all.html"
    umapTopPlotPath="reports/"+thisSampleID+"_"+thisSampleRef+"_UMAP_top.png"
    ht="<html><body><tt>sample ID: "+thisSampleID+" with reference "+thisSampleRef+"</tt><br>"
    ht=ht+"<a href='"+cnvPlotPath+"' target='_blank'><img align='Top' src='"+cnvPlotPath+"' width='50%' alt='CNV plot will appear here'></a>"
    ht=ht+"<a href='"+umapAllPlotlyPath+"' target='_blank'><img align='Top' src='"+umapAllPlotPath+"' width='50%' alt='UMAP plot will appear here'></a>"
    ht=ht+"</tt></table><body></html>"
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




class UI(object):
    """The CherryPy Web UI Webserver class defines entrypoints and
    function calls.
    """
    # global variables within the CherryPy Web UI
    cpgQueue = 0 # TODO use mutex instead
    umapQueue = 0
    cnvpQueue = 0
    cnv_lock = mp.Lock()
    umap_lock = mp.Lock()
    cpg_sem = threading.Semaphore()


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
            "cpgs": UI.cpgQueue,
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
                UI.cpgQueue = 0
            if queue_name == "umap":
                UI.umapQueue = 0
            if queue_name == "cnvp":
                UI.cnvpQueue = 0
            html += queue_name + " queue reset"
        return html

    @cherrypy.expose
    def status(self):
        positions = [pos.name for pos in minion_positions()]
        device_activity = {pos:real_device_activity(pos) for pos in positions}
        active_runs = {pos:active_run(pos) for pos in positions}
        sample_id = {pos:run_sample_id(pos) for pos in positions}
        url_live_device = {
            pos:url_for(UI.live_device_status, device_id=pos)
            for pos in positions
        }
        url_analysis = {
            pos:url_for(UI.AnalysisStatusLive, deviceString=pos)
            for pos in positions
        }
        url_device = {
            pos:url_for(UI.live_device_status, device_id=pos)
            for pos in positions
        }
        return render_template(
            "status.html",
            positions=positions,
            device_activity=device_activity,
            url_live_device=url_live_device,
            url_analysis=url_analysis,
            url_device=url_device,
            active_runs=active_runs,
            sample_id=sample_id,
            mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
        )

    @cherrypy.expose
    def start(
        self,
        device_id="",
        sample_id="",
        run_duration="",
        reference_id="",
        start_voltage="",
    ):
        start_now = bool(sample_id) and float(run_duration) >= 0.1
        if start_now:
            run_ids = start_run(
                device_id=device_id,
                sample_id=sample_id,
                run_duration=run_duration,
                start_voltage=start_voltage,
            )
            # TODO if not run_ids:
            write_reference_name(sample_id, reference_id)
            return render_template(
                "start.html",
                start_now=start_now,
                test=False,
                sample_id=sample_id,
                reference_id=reference_id,
                device_id=device_id,
                run_id=" / ".join(run_ids),
                mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
                run_info=run_information(device_id),
            )
        else:
            positions = [p.name for p in minion_positions()]
            idle = [p for p in positions if real_device_activity(p) == "idle"
                and flow_cell_id(p) != ""]
            flow_cell = {pos:flow_cell_id(pos) for pos in idle}
            return render_template(
                "start.html",
                start_now=start_now,
                test=False,
                idle=idle,
                flow_cell=flow_cell,
                references=reference_annotations(),
            )

    @cherrypy.expose
    def start_test(self, device_id=""):
        if device_id:
            sample_id = (date_time_string_now() + "_TestRun_"
                + flow_cell_id(device_id))
            run_ids = start_run(
                device_id=device_id,
                sample_id=sample_id,
                run_duration="0.1",
                start_voltage="-180",
            )
            # TODO if not run_ids:
            return render_template(
                "start.html",
                start_now=True,
                sample_id=sample_id,
                reference_id="TEST",
                device_id=device_id,
                run_id=" / ".join(run_ids),
                mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
                run_info=run_information(device_id),
            )
        else:
            positions = [p.name for p in minion_positions()]
            idle = [p for p in positions if real_device_activity(p) == "idle"
                and flow_cell_id(p) != ""]
            flow_cell = {pos:flow_cell_id(pos) for pos in idle}
            return render_template(
                "start.html",
                start_now=False,
                test=True,
                idle=idle,
                flow_cell=flow_cell,
                references=reference_annotations(),
            )

    @cherrypy.expose
    def stop_sequencing(self, device_id=""):
        protocol_id = stop_run(device_id)
        if protocol_id is None:
            return "No protocol running, nothing was stopped."
        else:
            return f"Protocol {protocol_id} stopped on {device_id}."

    @cherrypy.expose
    def list_runs(self):
        mounted_flow_cell_id = {}
        current_status = {}
        flow_cell_id = {}
        run_ids = {}
        device_names = []

        for minion in minion_positions():
            name = minion.name
            connection = minion.connect()
            device_names.append(name)
            mounted_flow_cell_id[name] = connection.device.get_flow_cell_info(
                ).flow_cell_id
            # READY, STARTING, sequencing/mux = PROCESSING, FINISHING;
            # Pause = PROCESSING
            current_status[name] = connection.acquisition.current_status()
            protocols = connection.protocol.list_protocol_runs()
            run_ids[name] = protocols.run_ids
            for run_id in run_ids[name]:
                run_info = connection.protocol.get_run_info(run_id=run_id)
                flow_cell_id[(name, run_id)] = run_info.flow_cell.flow_cell_id

        return render_template(
            "list_runs.html",
            device_names=device_names,
            host=CHERRYPY_HOST,
            mounted_flow_cell_id=mounted_flow_cell_id,
            current_status=current_status,
            flow_cell_id=flow_cell_id,
            run_ids=run_ids,
        )

    @cherrypy.expose
    def results(self):
        files = get_all_results()
        return render_template("results.html", files=files)

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
        if func == "cpg":
            return render_template(
                "analysis_cpg.html",
                start_time=date_time_string_now(),
                sample_name=samp,
                autorefresh="",
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
            args=(cnv_plt_data, UI.cnv_lock),
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
            args=(umap_data, UI.umap_lock),
            name="umap",
        )
        proc.start()
        proc.join()
        umap_data.read_from_disk()

        if close_up == "True":
            return umap_data.cu_plot_json

        return umap_data.plot_json

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
    def live_device_status(self,device_id=""):
        is_sequencing = real_device_activity(device_id) == "sequencing"
        sample_id = run_sample_id(device_id)
        state = run_state(device_id)
        yield_ = run_yield(device_id)
        status = None
        try:
            status = device_status(device_id)
            previous_activity = True
        except Exception as e:
            # TODO catch correct exception.
            print(e)
            sys.exit(1) # TODO del
            previous_activity = False
        return render_template(
            "device_status.html",
            device_id=device_id,
            status=status,
            flow_cell_id=flow_cell_id(device_id),
            is_sequencing=is_sequencing,
            sample_id=sample_id,
            yield_=yield_,
            state=state,
            previous_activity=previous_activity,
        )

    @cherrypy.expose
    def AnalysisStatusLive(self,deviceString=""):
        myString=""
        if deviceString:
            myString=livePage(deviceString)
        return myString

    @cherrypy.expose
    def cpgs(self, sample_name=""):
        """Generate a self-refreshing page to invoke methylation calling."""
        UI.cpg_sem.acquire()
        stats = methylation_caller(sample_name)
        UI.cpg_sem.release()
        return json.dumps(stats)

    @cherrypy.expose
    def launchAutoTerminator(self,sampleName="",deviceString=""):
        myString="ERROR"
        if sampleName and deviceString:
            myString=thisRunWatcherTerminator(deviceString,sampleName)
        return myString

def main():
    """Start CherryPy Webserver."""
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
    cherrypy.quickstart(UI(), "/", cherrypy_config)

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
