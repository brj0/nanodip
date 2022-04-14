"""
### CherryPy Web UI
The browser-based user interface is based on CherryPy, which contains an
intergrated web server and serves pages locally. Communication between the
service and browser typically generates static web pages that may or may not
contain automatic self refresh commands. In the case of self-refreshing pages,
the browser will re-request a given page with leads to re-execution of the
respective python functions. The main handles to these function are located in
the Web UI cell below.
"""

# start_external_modules
import cherrypy
import json
import pandas as pd
import psutil
import shutil
import socket
import os
import sys
import time
import threading
import multiprocessing as mp
# end_external_modules

# start_internal_modules
from config import (
    ANALYSIS_EXCLUSION_PATTERNS,
    BROWSER_FAVICON,
    CHERRYPY_HOST,
    CHERRYPY_PORT,
    DATA,
    DEBUG_MODE,
    ENDINGS,
    IMAGES,
    NANODIP_REPORTS,
    NEEDED_NUMBER_OF_BASES,
    VERBOSITY,
)
from utils import (
    convert_html_to_pdf,
    composite_path,
    date_time_string_now,
    reference_annotations,
    get_runs,
    predominant_barcode,
    render_template,
    read_reference,
    url_for,
    write_reference_name,
)
from api import (
    active_run,
    called_bases,
    device_status,
    flow_cell_id,
    get_all_results,
    minion_positions,
    methylation_caller,
    real_device_activity,
    run_information,
    run_sample_id,
    run_state,
    run_yield,
    set_bias_voltage,
    start_run,
    stop_run,
)
from data import (
    binary_reference_data_exists,
    ReferenceGenome,
)
from plots import (
    CNVData,
    UMAPData,
)
# end_internal_modules

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
    """List of Device, used to store auto-termination status for all
    devices.
    """
    def __init__(self):
        self.list = []

    def get(self, device_id):
        device = self.get_device(device_id)
        if device:
            return device
        else:
            device = Device(device_id)
            self.list.append(device)
            return device

    def pop(self, device_id):
        device = self.get_device(device_id)
        if device:
            self.list.remove(device)
            return device
        else:
            return None

    def get_device(self, device_id):
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

class UI(object):
    """The CherryPy Web UI Webserver class defines entrypoints and
    function calls.
    """
    # global variables within the CherryPy Web UI
    # TODO use Semaphore instead
    cnv_lock = mp.Lock()
    # TODO use Semaphore instead
    umap_lock = mp.Lock()
    cpg_sem = threading.Semaphore()
    devices = Devices()

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
            "cpgs": 1 - UI.cpg_sem._value,
            "cnvp": len([p for p in mp.active_children() if p.name == "cnv"]),
            "umap": len([p for p in mp.active_children() if p.name == "umap"]),
        }
        # Calculate urls to avoid hardcoding urls in html templates.
        return render_template(
            "index.html",
            sys_stat=sys_stat,
            url_cpgs=url_for(UI.reset_queue, queue_name=sys_stat["cpgs"]),
            url_cnvp=url_for(UI.reset_queue, queue_name=sys_stat["cnvp"]),
            url_umap=url_for(UI.reset_queue, queue_name=sys_stat["umap"]),
        )

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
    def restart(self):
        cherrypy.engine.restart()
        return render_template("restart.html")

    @cherrypy.expose
    def status(self):
        positions = [pos.name for pos in minion_positions()]
        device_activity = {}
        # Calculate urls to avoid hardcoding urls in html templates.
        url_analysis = {}
        url_device = {"none": UI.live_device_status.__name__}
        for pos in positions:
            device_activity[pos] = real_device_activity(pos)
            url_analysis[pos] = url_for(UI.live_plots, device_id=pos)
            url_device[pos] = url_for(UI.live_device_status, device_id=pos)
        return render_template(
            "status.html",
            positions=positions,
            url_analysis=url_analysis,
            url_device=url_device,
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
            # Delete termination status info of last round.
            UI.devices.pop(device_id)
            run_ids = start_run(
                device_id=device_id,
                sample_id=sample_id,
                run_duration=run_duration,
                start_voltage=start_voltage,
            )
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
        urls = {f:f"reports/{f}" for f in files}
        return render_template(
            "results.html",
            files=files,
            urls=urls,
        )

    @cherrypy.expose
    def analysis(self, func="", samp="", ref="", new="False"):
        if func == "":
            analysis_runs = [run for run in get_runs() if not any(
                pattern in run for pattern in ANALYSIS_EXCLUSION_PATTERNS)
            ]
            annotations = reference_annotations()
            # Calculate urls to avoid hardcoding urls in html templates.
            url_cnv = {}
            url_cnv_new = {}
            url_cpgs = {}
            url_pdf = {}
            url_umap = {}
            url_umap_new = {}
            for run in analysis_runs:
                url_cnv[run] = url_for(UI.analysis, func="cnv", samp=run)
                url_cnv_new[run] = url_for(
                    UI.analysis, func="cnv", samp=run, new=True,
                )
                url_cpgs[run] = url_for(UI.analysis, func="cpgs", samp=run)
                for a in annotations:
                    url_umap_new[(run,a)] = url_for(
                        UI.analysis, func="umap", samp=run, ref=a, new=True
                    )
                    url_umap[(run,a)] = url_for(
                        UI.analysis, func="umap", samp=run, ref=a
                    )
                    url_pdf[(run,a)] = url_for(
                        UI.make_pdf, samp=run, ref=a
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
            )
        if func == "cnv":
            genome = ReferenceGenome()
            genes = genome.genes.name.to_list()
            return render_template(
                "analysis_cnv.html",
                url_cnv=UI.cnv.__name__,
                sample_name=samp,
                genes=genes,
                new=new,
            )
        if func == "umap":
            return render_template(
                "analysis_umap.html",
                url_umap=UI.umap.__name__,
                sample_name=samp,
                reference_name=ref,
                new=new,
                first_use = not binary_reference_data_exists(),
            )
        if func == "cpgs":
            return render_template(
                "analysis_cpg.html",
                url_cpgs=UI.cpgs.__name__,
                start_time=date_time_string_now(),
                sample_name=samp,
                autorefresh="",
            )
        else:
            raise cherrypy.HTTPError(404, "URL not found")

    @cherrypy.expose
    def cnv(self, samp, genes="", new="False"):
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

        return cnv_plt_data.plot_cnv_and_genes([genes])

    @cherrypy.expose
    def umap(self, samp, ref, close_up="", new="False"):
        try:
            umap_data = UMAPData(samp, ref)
        except FileNotFoundError:
            raise cherrypy.HTTPError(405, "URL not allowed")

        def make_plot(plt_data, lock):
            """Plot function for multiprocessing."""
            lock.acquire()
            if not plt_data.files_on_disk() or new == "True":
                try:
                    plt_data.make_umap_plot()
                except ValueError:
                    raise cherrypy.HTTPError(405, "No data to plot.")
                    
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

    @cherrypy.expose
    def make_pdf(self, samp=None, ref=None):
        path = composite_path(NANODIP_REPORTS, samp, "cpgcount.txt")
        try:
            with open(path, "r") as f:
                overlap_cnt = f.read()
        except FileNotFoundError:
            raise cherrypy.HTTPError(
                405, 
                "CpG count file not found. Probably UMAP not completed."
            )
        path = composite_path(NANODIP_REPORTS, samp, ENDINGS["aligned_reads"])
        try:
            with open(path, "r") as f:
                read_numbers = f.read()
        except FileNotFoundError:
            raise cherrypy.HTTPError(
                405, 
                "Aligned reads file not found. Probably UMAP not completed."
            )

        cnv_path = composite_path(NANODIP_REPORTS, samp, ENDINGS["cnv_png"])
        umap_path = composite_path(
            NANODIP_REPORTS, samp, ref, ENDINGS["umap_top_png"],
        )
        pie_chart_path = composite_path(
            NANODIP_REPORTS, samp, ref, ENDINGS["pie"]
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
            pie_chart_path=pie_chart_path,
        )
        report_path = composite_path(
            NANODIP_REPORTS, samp, ref, ENDINGS["report"],
        )
        server_report_path = composite_path(
            "reports", samp, ref, ENDINGS["report"],
        )
        convert_html_to_pdf(html_report, report_path)
        raise cherrypy.HTTPRedirect(server_report_path)

    @cherrypy.expose
    def about(self):
        return render_template("about.html")

    @cherrypy.expose
    def live_device_status(self, device_id=""):
        UI.launch_auto_terminator(device_id)
        is_sequencing = real_device_activity(device_id) == "sequencing"
        sample_id = run_sample_id(device_id)
        state = run_state(device_id)
        yield_ = run_yield(device_id)
        status = None
        device = UI.devices.get(device_id)
        is_active = active_run(device_id) != "none" 
        try:
            status = device_status(device_id)
            previous_activity = True
        except Exception as e:
            # TODO catch correct exception.
            print("Name of the catched exception:", e)
            sys.exit(1)
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
            needed_mega_bases=NEEDED_NUMBER_OF_BASES // 1e6,
            url_auto_terminator = {"none": UI.set_auto_terminate.__name__},
            termination_type=device.termination_type,
            is_active=is_active,
        )

    @cherrypy.expose
    def set_auto_terminate(self, device_id="", termination_type=""):
        """Sets auto termination status."""
        device = UI.devices.get(device_id)
        if termination_type in ["terminated", "manually", "auto"]:
            device.termination_type = termination_type
            print(
                f"Auto terminate status of {device_id} "
                f"set to {device.termination_type}"
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
            called_bases(device.id) > NEEDED_NUMBER_OF_BASES
        ):
            stop_run(device.id)

    @cherrypy.expose
    def live_plots(self, device_id=""):
        """Generate a live preview of the data analysis with the current
        plots.
        """
        if not device_id:
            return ""
        # if there is a run that produces data, the run ID will exist
        sample_id = run_sample_id(device_id)
        reference = read_reference(sample_id)
        cnv_plt_path_png = composite_path(
            "reports", sample_id, ENDINGS["cnv_png"],
        )
        umap_plt_path_png = composite_path(
            "reports",
            sample_id, reference, ENDINGS["umap_all_png"],
        )
        umap_plt_path_html = composite_path(
            "reports",
            sample_id, reference, ENDINGS["umap_all_html"],
        )
        return render_template(
            "live_plots.html",
            sample_id=sample_id,
            reference=reference,
            cnv_plt_path_png=cnv_plt_path_png,
            umap_plt_path_png=umap_plt_path_png,
            umap_plt_path_html=umap_plt_path_html,
        )

    @cherrypy.expose
    def cpgs(self, sample_name=""):
        """Generate a self-refreshing page to invoke methylation calling."""
        UI.cpg_sem.acquire()
        stats = methylation_caller(sample_name)
        UI.cpg_sem.release()
        return json.dumps(stats)

    @cherrypy.expose
    def change_voltage(self, device_id="", voltage=""):      
        set_bias_voltage(device_id, voltage)
        return render_template(voltage=voltage)

def start_webserver():
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

