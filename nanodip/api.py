"""
### MinKNOW API Functions
Check https://github.com/nanoporetech/minknow_api for reference.

The following code requires a patched version of the MinKNOW API. Install it
from https://github.com/neuropathbasel/minknow_api.
"""

# start_external_modules
import argparse
import logging
import os
import shutil
import subprocess
import sys

from minknow_api.acquisition_pb2 import READY, STARTING, PROCESSING, FINISHING
from minknow_api.manager import Manager
from minknow_api.tools import protocols
import grpc
# end_external_modules

# start_internal_modules
from nanodip.config import (
    DATA,
    ENDING,
    F5C,
    MINIMAP2,
    NANODIP_OUTPUT,
    NEEDED_NUMBER_OF_BASES,
    READS_PER_FILE,
    REFERENCE_GENOME_FA,
    REFERENCE_GENOME_MMI,
    SAMTOOLS,
    THIS_HOST,
)
from nanodip.utils import (
    date_time_string_now,
    extract_referenced_cpgs,
    predominant_barcode,
    files_by_ending,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)

def minknow_manager():
    """Construct a manager using the host and port provided. This is
    used to connect to the MinKNOW service trough the MK API.

    minknow_api.manager.Manager:  a wrapper around MinKNOW's Manager
        gRPC API with utilities for querying sequencing positions and
        offline basecalling tools.
    """
    return Manager(host=THIS_HOST, port=9501, use_tls=False)

def minion_positions():
    """Return MinION devices that are currenty connected to the system."""
    manager = minknow_manager()
    # Find a list of currently available sequencing positions.
    positions = manager.flow_cell_positions()
    # User could call {posisions.connect()} here to connect to the
    # running MinKNOW instance.
    return positions

def connection_from_device_id(device_id):
    """Returns minion position."""
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
    cell_id = "no_flow_cell"
    positions = minion_positions()
    for p in positions:
        if device_id in p.name:
            connection = p.connect()
            cell_id = connection.device.get_flow_cell_info().flow_cell_id
    return cell_id

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
    # 'Read until' must have a reference and a filter type, if enabled:
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
        "--",  # Required for so-called extra-arguments.
        "--start_bias_voltage", start_voltage,
    ]

    # Parse arguments to be passed to started protocols:
    args = parse_args(args_list)

    # # Specify --verbose on the command line to get extra details.
    # if args.verbose:
        # logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

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
    protocol = connection.protocol.list_protocol_runs()
    protocol_id = protocol.run_ids[-1]
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

def device_activity(device_id):
    """Returns device activity. Virtual test runs will be recognized as
    active.
    """
    connection = connection_from_device_id(device_id)
    status = connection.acquisition.current_status().status
    device_activity = {
        STARTING: "sequencing",
        PROCESSING: "sequencing",
        FINISHING: "sequencing",
        READY: "idle",
    }
    return device_activity.get(status, "")

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

def number_of_called_bases(device_id):
    """Returns number of called bases."""
    connection = connection_from_device_id(device_id)
    # Check if device is working.
    if connection.acquisition.current_status().status == READY:
        return 0
    acquisition = connection.acquisition.get_acquisition_info()
    num_bases = acquisition.yield_summary.estimated_selected_bases
    return num_bases

def device_status(device_id):
    """MinKNOW status for device {device_id}."""
    connection = connection_from_device_id(device_id)
    current_bases = number_of_called_bases(device_id)
    needed_mb = round(NEEDED_NUMBER_OF_BASES // 1e6, 2)
    current_mb = round(current_bases / 1e6, 2)
    progress = round(100*current_mb/needed_mb, 1)
    status = {
        "Real device activity": real_device_activity(device_id),
        "Active run": active_run(device_id),
        "Progress": f"{progress}% ({current_mb} MB / {needed_mb} MB)",
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
    # Progress is only needed if device is sequencing.
    if active_run(device_id) == "none":
        status.pop("Progress")
    return status

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
    """Get sample ID from MinKNOW by device, only available after data
    acquisition has been initiated by MinKNOW.
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

def set_bias_voltage(device_id, voltage):
    """Change MinKnow bias voltage."""
    connection = connection_from_device_id(device_id)
    previous_voltage = connection.device.get_bias_voltage().bias_voltage
    connection.device.set_bias_voltage(bias_voltage=float(voltage))

def single_file_methylation_caller(analysis_dir):
    """Invokes f5c methylation caller on a single fast5/fastq file and
    calculates methylation frequencies and CpG overlaps.

    Args:
        analysis_dir: Directory containing fast5 and fastq files to be
            analyzed. The basename of these files (corresponding to the
            run id) must be equal to the name of analysis_dir.
            The resulting files of the methylation calling process will
            be saved in analysis_dir.
    """
    file_name = os.path.basename(analysis_dir)
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
        "-o", base_path + ENDING["readsort_bam"],
    ]
    # Make bam index for samtools.
    bam_index = [
        SAMTOOLS, "index",
        base_path + ENDING["readsort_bam"],
    ]
    # Methylation caller.
    methyl_calling = [
        F5C, "call-methylation",
        #"--disable-cuda=yes",   # For debugging on CPU only.
        "-B2000000", "-K400",    # Set B to 2 megabases (GPU) and 0.4 kreads
        "-b", base_path + ENDING["readsort_bam"],
        "-g", REFERENCE_GENOME_FA,
        "-r", base_path + ".fastq",
        ">", base_path + ENDING["result_tsv"],
    ]
    # Calculate methylation frequencies.
    methyl_frequency = [
        F5C, "meth-freq",
        "-c", "2.5",
        "-s",
        "-i", base_path + ENDING["result_tsv"],
        ">",
        base_path + ENDING["freq_tsv"],
    ]
    commands = [
        f5c_index,
        seq_align,
        bam_index,
        methyl_calling,
        methyl_frequency,
    ]

    def log_subprocess_output(pipe):
        """Logs subprocess (to file) and also sends output to stdout."""
        stdout_data, stderr_data = pipe.communicate()
        for line in stdout_data.decode().split("\n"):
            logger.info(line)
            print(line)
        for line in stderr_data.decode().split("\n"):
            logger.info(line)
            print(line)

    for cmd in commands:
        cmd_str = " ".join(cmd)
        process = subprocess.Popen(
            cmd_str,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        log_subprocess_output(process)
        exit_code = process.wait()
        if exit_code != 0:
            logger.error("Error occured on subprocess '%s'", cmd[0])
    # Calculate CpG overlap.
    extract_referenced_cpgs(
        base_path + ENDING["freq_tsv"],
        base_path + ENDING["methoverl_tsv"],
        base_path + ENDING["methoverlcnt_txt"],
    )
    print(f"Methylation calling on {file_name} done.")

def remove_dirs_with_wrong_barcode(sample_name, barcode):
    """Removes directories with wrong barcode due to change of
    the predominant barcode after the first reads.
    """
    try:
        output_dirs = [
            os.path.join(NANODIP_OUTPUT, sample_name, dir_nm) for dir_nm in
            os.listdir(os.path.join(NANODIP_OUTPUT, sample_name))
        ]
    except FileNotFoundError:
        return
    wrong_barcode_dirs = [
        f for f in output_dirs if barcode not in f
    ]
    for dir_path in wrong_barcode_dirs:
        shutil.rmtree(dir_path)

def methylation_calling_done(analysis_dir):
    """Checks if methylation calling is done."""
    endings = [
        ENDING[e]
        for e in [
            "methoverlcnt_txt",
            "methoverl_tsv",
            "result_tsv",
            "freq_tsv",
            "readsort_bai",
            "readsort_bam",
        ]
    ]
    ending_present = [
        any(file_.endswith(ending) for file_ in os.listdir(analysis_dir))
        for ending in endings
    ]
    return all(ending_present)

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
    remove_dirs_with_wrong_barcode(sample_name, barcode)
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

    # Collect all passed fast5/fastq pairs
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
        if methylation_calling_done(analysis_dir):
            prev_called.append(file_name)
        else:
            not_called.append(
                [analysis_dir, file_name]
            )
    for directory, file_name in not_called:
        single_file_methylation_caller(directory)
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
