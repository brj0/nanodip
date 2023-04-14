"""
### CherryPy Web UI
The browser-based user interface is based on CherryPy, which contains an
integrated web server and serves pages locally. Communication between the
service and browser typically generates static web pages that may or may not
contain automatic self refresh commands. In the case of self-refreshing pages,
the browser will re-request a given page with leads to re-execution of the
respective python functions. The main handles to these function are located in
the Web UI cell below.
"""

# start_external_modules
import json
import logging
import os
import re
import shutil
import socket
import threading

import cherrypy
import grpc
import psutil
# end_external_modules

# start_internal_modules
from nanodip.config import (
    ANALYSIS_EXCLUSION_PATTERNS,
    CHERRYPY_HOST,
    CHERRYPY_PORT,
    DATA,
    DEBUG_MODE,
    ENDING,
    EPIDIP_UMAP_COORDINATE_FILES,
    NANODIP_REPORTS,
    NEEDED_NUMBER_OF_BASES,
)
from nanodip.utils import (
    composite_path,
    convert_html_to_pdf,
    date_time_string_now,
    get_all_results,
    get_runs,
    predominant_barcode,
    reference_annotations,
    render_template,
    url_for,
)
from nanodip.api import (
    active_run,
    device_status,
    flow_cell_id,
    methylation_caller,
    minion_positions,
    number_of_called_bases,
    real_device_activity,
    run_information,
    run_sample_id,
    run_state,
    run_yield,
    set_bias_voltage,
    start_run,
    stop_run,
)
from nanodip.data import (
    Genome,
    Reference,
    Sample,
    binary_reference_data_exists,
)
from nanodip.plots import (
    CNVData,
    UMAPData,
)
from nanodip.epidip import(
    download_epidip_data,
    top_variable_cpgs,
)
from nanodip.classifiers import (
    fit_and_evaluate_classifiers,
)
# end_internal_modules

# Define logger
logger = logging.getLogger(__name__)

class Device:
    """Container used to store auto-termination status for a single
    device.
    """
    def __init__(self, device_id):
        self.id = device_id
        self.termination_type = "manually"

    def __repr__(self):
        return f"(id={self.id}, termination_type={self.termination_type})"

class Devices:
    """List of Device objects, used to store auto-termination status for
    all devices.
    """
    def __init__(self):
        self.list = []

    def get(self, device_id):
        """Appends device {device_id} if necessary and returns it."""
        device = self.get_device(device_id)
        if device:
            return device
        device = Device(device_id)
        self.list.append(device)
        return device

    def pop(self, device_id):
        """Removes and returns device."""
        device = self.get_device(device_id)
        if device:
            self.list.remove(device)
            return device
        return None

    def get_device(self, device_id):
        """Returns device with given id. Returns false if id is not found."""
        return next(
            (device for device in self.list if device_id == device.id),
            False,
        )

    def __iter__(self):
        return iter(self.list)

    def __contains__(self, other):
        return other in [d.id for d in self]

    def __repr__(self):
        out = "["
        out += ", ".join([str(d) for d in self.list])
        out += "]"
        return out

class ActivePlots:
    """Container to keep the most recent plots in memory to speed up
    computation for CNV gene plot.
    """
    def __init__(self):
        self.plot = {}
    def _delete_old_plots(self):
        """Limits plot dictionary to 5 elements and deletes oldest plots
        first.
        """
        while len(self.plot) > 5:
            # Since dict is insertion ordered this is the oldest entry.
            old_id = list(self.plot.keys())[0]
            self.plot.pop(old_id)
    def set(self, id_, active_plot):
        """Saves plot data with given identifier."""
        self.plot[id_] = active_plot
        self._delete_old_plots()
    def get(self, webpage_id):
        """Returns CNVData with given id."""
        active_plot =self.plot.get(webpage_id, False)
        if active_plot:
            return active_plot
        return None
    def __repr__(self):
        """Enables print and str for debugging purposes."""
        return str(self.plot)

class MultiSemaphore:
    """Two-level semaphore object allowing one process using light and
    one process using heavy utilization of the system simultaneously.
    """

    def __init__(self):
        self.light = threading.Semaphore()
        self.heavy = threading.Semaphore()
        self.proc_cnt = {"cnv": 0, "umap": 0, "cpg": 0, "clf": 0}
        self.workload = {
            "cnv": "light",
            "umap": "heavy",
            "cpg": "heavy",
            "clf": "heavy",
        }

    def acquire(self, type_):
        """Acquire lock."""
        if self.workload[type_] == "light":
            self.light.acquire()
            self.proc_cnt[type_] += 1
        elif self.workload[type_] == "heavy":
            self.heavy.acquire()
            self.proc_cnt[type_] += 1
        else:
            raise ValueError("Invalid type.")

    def release(self, type_):
        """Release lock."""
        if self.workload[type_] == "light":
            self.proc_cnt[type_] -= 1
            self.light.release()
        elif self.workload[type_] == "heavy":
            self.proc_cnt[type_] -= 1
            self.heavy.release()
        else:
            raise ValueError("Invalid type.")

    def __str__(self):
        """Prints overview of object for debugging purposes."""
        return "Number of running processes:\n" + str(self.proc_cnt)


class UI:
    """User interface implemented as CherryPy webserver."""
    # global variables within the CherryPy Web UI
    devices = Devices()
    active_plots = ActivePlots()
    sem = MultiSemaphore()

    # Methylation calling, UMAP and and non-supervised classifiers all
    # heavily use system resources, and using them at the same time
    # leads to crash. Therefore a shared semaphore is used.
    all_sem = threading.Semaphore()

    @cherrypy.expose
    def index(self, refresh="False"):
        """Start page."""
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
            "cpgs": UI.sem.proc_cnt["cpg"],
            "cnvp": UI.sem.proc_cnt["cnv"],
            "umap": UI.sem.proc_cnt["umap"],
            "clf": UI.sem.proc_cnt["clf"],
        }
        # Calculate URL to avoid hard coding URLs in HTML templates.
        return render_template(
            "index.html",
            sys_stat=sys_stat,
            url_restart=url_for(UI.restart),
            url_refresh=url_for(UI.index, refresh="True"),
            refresh=refresh,
        )

    @cherrypy.expose
    def restart(self):
        """Restart CherryPy."""
        cherrypy.engine.restart()
        return render_template("restart.html")

    @cherrypy.expose
    def status(self):
        """List all devices with run-status and show UMAP/CMV preview if
        available.
        """
        device_ids = [pos.name for pos in minion_positions()]
        # Calculate URLs to avoid hard coding URLs in HTML templates.
        return render_template(
            "status.html",
            device_ids=device_ids,
            url_live_device_status=url_for(UI.status_device),
            url_live_plots=url_for(UI.status_plots),
        )

    @cherrypy.expose
    def status_device(self, device_id):
        """Lists run-status of device {device_id} and prints button for
        terminating run.
        """
        UI.launch_auto_terminator(device_id)
        status = None
        device = UI.devices.get(device_id)
        is_active = active_run(device_id) != "none"
        try:
            status = device_status(device_id)
            previous_activity = True
        except grpc._channel._InactiveRpcError:
            previous_activity = False
        return render_template(
            "status_device.html",
            device_id=device_id,
            status=status,
            flow_cell_id=flow_cell_id(device_id),
            sample_id=run_sample_id(device_id),
            yield_=run_yield(device_id),
            state=run_state(device_id),
            previous_activity=previous_activity,
            needed_mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
            url_auto_terminator = url_for(UI.set_auto_terminate),
            termination_type=device.termination_type,
            is_active=is_active,
        )

    @cherrypy.expose
    def status_plots(self, device_id=""):
        """Generate a live preview of the data analysis with the current
        plots.
        """
        if not device_id:
            raise cherrypy.HTTPError(404, "URL not found")
        # If there is a run that produces data, the run ID will exist.
        sample_id = run_sample_id(device_id)

        # Chose reference with latest png-UMAP file.
        umap_png_files = [
            f for f in os.listdir(NANODIP_REPORTS)
            if f.endswith(ENDING["umapall_png"]) and sample_id in f
        ]
        if not umap_png_files:
            # Dummy reference name if no UMAP files are found (happens if
            # device is not sequencing and thus no sample_id can be found.
            reference = "none"
        else:
            latest_umap_png = max(
                [os.path.join(NANODIP_REPORTS, f) for f in umap_png_files],
                key=os.path.getmtime,
            )
            reference = re.search(
                sample_id +  "_(.*?)_" + ENDING["umapall_png"],
                latest_umap_png
            ).group(1)

        cnv_plt_path_png = composite_path(
            "reports", sample_id, ENDING["cnv_png"],
        )
        cnv_plt_path_html = composite_path(
            "reports", sample_id, ENDING["cnv_html"],
        )
        umap_plt_path_png = composite_path(
            "reports", sample_id, reference, ENDING["umapall_png"],
        )
        umap_plt_path_html = composite_path(
            "reports", sample_id, reference, ENDING["umapall_html"],
        )
        return render_template(
            "status_plots.html",
            sample_id=sample_id,
            reference=reference,
            cnv_plt_path_png=cnv_plt_path_png,
            cnv_plt_path_html=cnv_plt_path_html,
            umap_plt_path_png=umap_plt_path_png,
            umap_plt_path_html=umap_plt_path_html,
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
        """Start sequencing run."""
        start_now = bool(sample_id) and float(run_duration) >= 0.1
        if start_now:
            # Delete termination status info of last round.
            UI.devices.pop(device_id)
            run_ids = start_run(
                device_id=device_id,
                sample_id=sample_id,
                run_duration=run_duration,
                start_voltage=start_voltage,
            )
            return render_template(
                "start.html",
                url_action=url_for(UI.start),
                start_now=start_now,
                test=False,
                sample_id=sample_id,
                reference_id=reference_id,
                device_id=device_id,
                run_id=" / ".join(run_ids),
                run_info=run_information(device_id),
            )
        positions = [p.name for p in minion_positions()]
        idle = [
            p for p in positions if real_device_activity(p) == "idle"
            and flow_cell_id(p) != ""
        ]
        flow_cell = {pos:flow_cell_id(pos) for pos in idle}
        return render_template(
            "start.html",
            url_action=url_for(UI.start),
            start_now=start_now,
            test=False,
            idle=idle,
            flow_cell=flow_cell,
            references=reference_annotations(),
        )

    @cherrypy.expose
    def start_test(self, device_id=""):
        """Start sequencing test run."""
        if device_id:
            sample_id = (date_time_string_now() + "_TestRun_"
                + flow_cell_id(device_id))
            run_ids = start_run(
                device_id=device_id,
                sample_id=sample_id,
                run_duration="0.1",
                start_voltage="-180",
            )
            return render_template(
                "start.html",
                start_now=True,
                sample_id=sample_id,
                reference_id="TEST",
                device_id=device_id,
                run_id=" / ".join(run_ids),
                run_info=run_information(device_id),
            )
        positions = [p.name for p in minion_positions()]
        idle = [p for p in positions if real_device_activity(p) == "idle"
            and flow_cell_id(p) != ""]
        flow_cell = {pos:flow_cell_id(pos) for pos in idle}
        return render_template(
            "start.html",
            url_action=url_for(UI.start_test),
            start_now=False,
            test=True,
            idle=idle,
            flow_cell=flow_cell,
            references=reference_annotations(),
        )

    @cherrypy.expose
    def stop_sequencing(self, device_id=""):
        """Can be used to manually stop a run."""
        protocol_id = stop_run(device_id)
        if protocol_id is None:
            return "No protocol running, nothing was stopped."
        return f"Protocol {protocol_id} stopped on {device_id}."

    @cherrypy.expose
    def list_runs(self):
        """Lists running and buffered sequencing runs."""
        mounted_flow_cell_id = {}
        current_status = {}
        flow_cell = {}
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
                flow_cell[(name, run_id)] = run_info.flow_cell.flow_cell_id

        return render_template(
            "list_runs.html",
            device_names=device_names,
            host=CHERRYPY_HOST,
            mounted_flow_cell_id=mounted_flow_cell_id,
            current_status=current_status,
            flow_cell=flow_cell,
            run_ids=run_ids,
        )

    @cherrypy.expose
    def results(self):
        """Lists all files in NANODIP_REPORTS."""
        files = get_all_results()
        urls = {f:f"reports/{f}" for f in files}
        return render_template(
            "results.html",
            files=files,
            urls=urls,
        )

    @cherrypy.expose
    def analysis(
        self,
        func="",
        sample_name="",
        reference_name="",
        new="False",
    ):
        """Creates an overview of all samples and provides the analysis
        tools (CpG methylation calling, CNV plot and UMAP plot).
        """
        if func == "":
            analysis_runs = [
                run for run in get_runs() if not any(
                    pattern in run for pattern in ANALYSIS_EXCLUSION_PATTERNS
                )
            ]
            annotations = reference_annotations()
            # Calculate URLs to avoid hard coding URLs in HTML templates.
            url_cnv = {}
            url_cnv_new = {}
            url_cpgs = {}
            url_pdf = {}
            url_umap = {}
            url_umap_new = {}
            url_clf = {}
            for run in analysis_runs:
                url_cnv[run] = url_for(
                    UI.analysis, func="cnv", sample_name=run,
                )
                url_cnv_new[run] = url_for(
                    UI.analysis, func="cnv", sample_name=run, new=True,
                )
                url_cpgs[run] = url_for(
                    UI.analysis, func="cpgs", sample_name=run,
                )
                for annotation in annotations:
                    url_umap_new[(run, annotation)] = url_for(
                        UI.analysis,
                        func="umap",
                        sample_name=run,
                        reference_name=annotation,
                        new=True,
                    )
                    url_umap[(run, annotation)] = url_for(
                        UI.analysis,
                        func="umap",
                        sample_name=run,
                        reference_name=annotation,
                    )
                    url_pdf[(run, annotation)] = url_for(
                        UI.make_pdf,
                        sample_name=run,
                        reference_name=annotation,
                    )
                    url_clf[(run, annotation)] = url_for(
                        UI.classifiers,
                        sample_name=run,
                        reference_name=annotation,
                        start="True",
                    )
            return render_template(
                "analysis_start.html",
                analysis_runs=analysis_runs,
                annotations=annotations,
                url_cnv=url_cnv,
                url_cnv_new=url_cnv_new,
                url_cpgs=url_cpgs,
                url_pdf=url_pdf,
                url_umap=url_umap,
                url_umap_new=url_umap_new,
                url_clf=url_clf,
            )
        if func == "cnv":
            genome = Genome()
            genes = genome.genes.name.to_list()
            return render_template(
                "analysis_cnv.html",
                url_cnv=url_for(UI.cnv),
                sample_name=sample_name,
                genes=genes,
                new=new,
            )
        if func == "umap":
            return render_template(
                "analysis_umap.html",
                url_umap=url_for(UI.umap_plt),
                sample_name=sample_name,
                reference_name=reference_name,
                new=new,
                first_use = not binary_reference_data_exists(),
            )
        if func == "cpgs":
            return render_template(
                "analysis_cpg.html",
                url_cpgs=url_for(UI.cpgs),
                start_time=date_time_string_now(),
                sample_name=sample_name,
            )
        raise cherrypy.HTTPError(404, "URL not found")

    @cherrypy.expose
    def cnv(self, sample_name, genes="", new="False", browser_tab_id=""):
        """Creates CNV plot and returns it as JSON.
            Args:
                sample_name: Name of sample to analyze.
                genes: List of genes that should be plotted as well.
                new: If true CNV is constructed de novo, otherwise data is
                    read from disk if possible.
                browser_tab_id: Unique identifier for each browser tab, used
                    to send data still in memory to correct front-end.
        """
        if genes and browser_tab_id:
            cnv_data = UI.active_plots.get(browser_tab_id)
            return cnv_data.plot_cnv_and_genes([genes])
        UI.sem.acquire("cnv")
        try:
            cnv_data = CNVData(sample_name)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")

        if not cnv_data.files_on_disk() or new == "True":
            try:
                cnv_data.make_cnv_plot()
            except ValueError:
                UI.sem.release("cnv")
                raise cherrypy.HTTPError(405, "No data to plot.")
        else:
            cnv_data.read_from_disk()
        UI.sem.release("cnv")
        UI.active_plots.set(browser_tab_id, cnv_data)
        return cnv_data.plot_cnv_and_genes()

    @cherrypy.expose
    def umap_plt(self, reference_name, sample_name, new="False", ntop=""):
        """Creates UMAP plot and returns it as JSON."""
        UI.sem.acquire("umap")
        try:
            umap_data = UMAPData.from_names(reference_name, sample_name)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")
        if not umap_data.files_on_disk() or new == "True":
            try:
                umap_data.make_umap_plot()
            except ValueError:
                UI.sem.release("umap")
                raise cherrypy.HTTPError(405, "No data to plot.")
        else:
            umap_data.read_from_disk()
        UI.sem.release("umap")
        return json.dumps({
            "all": umap_data.plot_json,
            "close_up": umap_data.cu_plot_json
        })

    @cherrypy.expose
    def umap_plt_epidip(
        self,
        reference_name,
        ntop,
        sample_name="",
        new="",
    ):
        """Creates UMAP plot and returns it as JSON."""
        try:
            nr_top_cpgs = int(ntop)
        except ValueError:
            raise cherrypy.HTTPError(405, "'ntop' must be an integer")
        UI.sem.acquire("umap")
        try:
            reference = Reference(reference_name)
            top_cpgs = top_variable_cpgs(reference, nr_top_cpgs)
            sample = Sample.from_cpgs(top_cpgs)
            umap_data = UMAPData(reference, sample)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")
        try:
            umap_data.make_umap_plot()
        except ValueError:
            UI.sem.release("umap")
            raise cherrypy.HTTPError(405, "No data to plot.")
        else:
            umap_data.read_from_disk()
        UI.sem.release("umap")
        return json.dumps({
            "all": umap_data.plot_json,
        })

    @cherrypy.expose
    def epidip_umap(
        self,
        reference_name="",
        ntop="",
        sample_name="",
    ):
        """Page for invoking UMAP for EpiDiP reference data."""
        if bool(reference_name) and bool(ntop):
            return render_template(
                "analysis_umap.html",
                url_umap=url_for(UI.umap_plt_epidip),
                sample_name=sample_name,
                reference_name=reference_name,
                ntop=ntop,
                first_use = not binary_reference_data_exists(),
            )
        return render_template(
            "epidip_umap.html",
            url_action=url_for(UI.epidip_umap),
            references=reference_annotations(),
        )

    @cherrypy.expose
    def classifiers(
        self,
        sample_name,
        reference_name,
        start="False",
    ):
        """Start non supervised classifiers."""
        return render_template(
            "classifiers.html",
            url_start_and_refresh=url_for(UI.classifiers_start_and_refresh),
            sample_name=sample_name,
            reference_name=reference_name,
            start=start,
        )

    @cherrypy.expose
    def classifiers_start_and_refresh(
        self,
        sample_name,
        reference_name,
        start="False",
    ):
        """Used to start training/evaluating classifiers and refresh
        evaluation output.
        """
        file_path = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["clf_txt"],
        )
        if start == "True":
            with open(file_path, "w") as f:
                f.write("Classification started. Results will appear here.")
            UI.sem.acquire("clf")
            fit_and_evaluate_classifiers(sample_name, reference_name)
            UI.sem.release("clf")
        try:
            with open(file_path, "r") as f:
                clf_results = f.read()
            return clf_results
        except FileNotFoundError:
            return "No Data to plot."

    @cherrypy.expose
    def make_pdf(self, sample_name=None, reference_name=None):
        """Generates PDF report."""
        path_cgp = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["cpg_cnt"]
        )
        path_reads = composite_path(
            NANODIP_REPORTS, sample_name, ENDING["alignedreads_txt"]
        )
        try:
            with open(path_cgp, "r") as f:
                overlap_cnt = f.read()
            with open(path_reads, "r") as f:
                read_numbers = f.read()
        except FileNotFoundError:
            overlap_cnt = 0
            read_numbers = 0
        cnv_path = composite_path(
            NANODIP_REPORTS, sample_name, ENDING["cnv_png"]
        )
        umap_path = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["umaptop_png"],
        )
        pie_chart_path = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["pie_png"]
        )
        logo_path = cherrypy.config["/favicon.ico"]["tools.staticfile.filename"]
        html_report = render_template(
            "pdf_report.html",
            sample_name=sample_name,
            logo_path=logo_path,
            sys_name=socket.gethostname(),
            date=date_time_string_now(),
            barcode=predominant_barcode(sample_name),
            reads=read_numbers,
            cpg_overlap_cnt=overlap_cnt,
            reference=reference_name,
            cnv_path=cnv_path,
            umap_path=umap_path,
            pie_chart_path=pie_chart_path,
        )
        report_path = composite_path(
            NANODIP_REPORTS, sample_name, reference_name, ENDING["report_pdf"],
        )
        server_report_path = composite_path(
            "reports", sample_name, reference_name, ENDING["report_pdf"],
        )
        convert_html_to_pdf(html_report, report_path)
        raise cherrypy.HTTPRedirect(server_report_path)

    @cherrypy.expose
    def set_auto_terminate(self, device_id="", termination_type=""):
        """Sets auto termination status."""
        device = UI.devices.get(device_id)
        if termination_type in ["terminated", "manually", "auto"]:
            device.termination_type = termination_type
            logger.info(
                "Auto terminate status of %s set to %s",
                device_id,
                device.termination_type,
            )
        else:
            raise cherrypy.HTTPError(
                404, "Invalid termination type: '{termination_type}'"
            )
        if termination_type == "terminated":
            stop_run(device.id)

    def launch_auto_terminator(device_id=""):
        """Terminates the current run if auto terminator is set."""
        device = UI.devices.get(device_id)
        if (
            device.termination_type == "auto" and
            number_of_called_bases(device.id) > NEEDED_NUMBER_OF_BASES
        ):
            stop_run(device.id)

    @cherrypy.expose
    def cpgs(self, sample_name=""):
        """Invokes methylation calling and returns statistics."""
        UI.sem.acquire("cpg")
        stats = methylation_caller(sample_name)
        UI.sem.release("cpg")
        return json.dumps(stats)

    @cherrypy.expose
    def change_voltage(self, device_id="", voltage=""):
        """Change bias voltage."""
        set_bias_voltage(device_id, voltage)
        return render_template(
            "change_voltage.html",
            voltage=voltage,
        )

    @cherrypy.expose
    def epidip_report(
        self,
        sentrix_id=None,
        reference_id=None,
        reference_umap=None,
    ):
        """Create report for reference case {sentrix_id} with precalculated
        umap coordinates.
        """
        if sentrix_id and reference_id and reference_umap:
            download_epidip_data(sentrix_id, reference_umap)
            umap_data = UMAPData.from_names(reference_id, sentrix_id)
            umap_data.read_precalculated_umap_matrix(reference_umap)
            umap_data.draw_pie_chart()
            umap_data.draw_scatter_plot()
            umap_data.draw_cu_scatter_plot()
            umap_data.save_to_disk()
            UI.make_pdf(
                self, sample_name=sentrix_id, reference_name=reference_id
            )
        else:
            return render_template(
                "epidip_report.html",
                url_action=url_for(UI.epidip_report),
                reference_umap=reference_umap,
                epidip_umaps=EPIDIP_UMAP_COORDINATE_FILES,
                references=reference_annotations(),
            )

    @cherrypy.expose
    def about(self):
        """About Page."""
        return render_template("about.html")


def start_webserver():
    """Start CherryPy Webserver."""
    if DEBUG_MODE:
        #Set access logging
        cherrypy.log.screen = True
        cherrypy.config.update({'log.screen': True})
    else:
        #Set access logging
        cherrypy.log.screen = False
        cherrypy.config.update({'log.screen': False})
        cherrypy.config.update({ "environment": "embedded" })
    cherrypy.log.access_file = composite_path(NANODIP_REPORTS, "cherrypy.log")
    cherrypy.log.error_file = composite_path(NANODIP_REPORTS, "cherrypy.log")

    print(f"NanoDiP server running at http://{CHERRYPY_HOST}:{CHERRYPY_PORT}")

    app_dir = os.path.dirname(__file__)
    cherrypy_config = {
        '/favicon.ico': {
            'tools.staticfile.on': True,
            'tools.staticfile.filename': (
                os.path.join(app_dir, "static/img/logo.png")
            ),
        },
        '/reports': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': NANODIP_REPORTS,
        },
        '/static': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': os.path.join(app_dir, "static"),
        },
    }
    cherrypy.quickstart(UI(), "/", cherrypy_config)
