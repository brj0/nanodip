"""
## Utils

Contains general utility functions.
"""


# start_external_modules
import os
import warnings

import colorsys
import datetime
import hashlib
import inspect
import jinja2
import pandas as pd
import xhtml2pdf.pisa
# end_external_modules

# start_internal_modules
from config import (
    ANNOTATIONS,
    ANNOTATIONS_ABBREVIATIONS_BASEL,
    ANNOTATIONS_ABBREVIATIONS_TCGA,
    BARCODE_NAMES,
    BETA_VALUES,
    BROWSER_FAVICON,
    CHROMOSOMES,
    DATA,
    ENDING,
    EXCLUDED_FROM_ANALYSIS,
    F5C,
    GENES,
    GENES_RAW,
    ILLUMINA_CG_MAP,
    MINIMAP2,
    NANODIP_REPORTS,
    REFERENCE_GENOME_FA,
    REFERENCE_GENOME_MMI,
    RELEVANT_GENES,
    RESULT_ENDING,
    SAMTOOLS,
)
# end_internal_modules

def sanity_check():
    """Checks if reference data is available and generates a warning
    otherwise.
    """
    requested_files = [
        BETA_VALUES,
        ANNOTATIONS_ABBREVIATIONS_BASEL,
        ANNOTATIONS_ABBREVIATIONS_TCGA,
        ILLUMINA_CG_MAP,
        CHROMOSOMES,
        REFERENCE_GENOME_FA,
        REFERENCE_GENOME_MMI,
        GENES_RAW,
        GENES,
        RELEVANT_GENES,
        BROWSER_FAVICON,
        F5C,
        MINIMAP2,
        SAMTOOLS,
    ]
    for f in requested_files:
        if not os.path.exists(f):
            warnings.warn(
                f"File '{f}' not found.\nFunctionality may be restricted.",
                RuntimeWarning,
            )

def extract_referenced_cpgs(sample_methylation,
                            output_overlap,
                            output_overlap_cnt):
    """Extract Illumina CpG sites including methylation status from sample.
    Sex chromosomes are removed.

    Args:
        sample_methylation: methylation file of sample
        output_overlap: file path of CpG overlap
        output_overlap_cnt: file path of CpG overlap count
    """
    reference_cpgs = pd.read_csv(
        ILLUMINA_CG_MAP,
        delimiter="\t",
        names=["ilmnid", "chromosome", "strand", "start"],
    )
    sample_cpgs = pd.read_csv(
        sample_methylation,
        delimiter="\t",
    )
    cpgs = pd.merge(sample_cpgs, reference_cpgs, on=["chromosome", "start"])
    # Extract singelton CpG's
    cpgs = cpgs.loc[cpgs["num_cpgs_in_group"] == 1]
    # Remove duplicates and sex chromosomes
    cpgs = cpgs.loc[
       (~cpgs["chromosome"].isin(["chrX", "chrY"]))
       & (~cpgs["ilmnid"].duplicated())
    ]
    cpgs["is_methylated"] = 0
    cpgs.loc[cpgs["methylated_frequency"] > 0.5, "is_methylated"] = 1
    # Write overlap Data Frame
    cpgs[["ilmnid", "is_methylated"]].to_csv(
        output_overlap, header=False, index=False, sep="\t",
    )
    # Write number of CpG's
    with open(output_overlap_cnt, "w") as f:
        f.write(f"{len(cpgs)}")

def render_template(template_name, **context):
    """Renders jinja2 templates to HTML."""
    loader = jinja2.FileSystemLoader("templates")
    template = jinja2.Environment(
        loader=loader).get_template(template_name)
    return template.render(context)

def url_for(url_func, **args):
    """Transforms a cherrypy function together with argument list to
    url string.
    Example:
        url_for(CherryPyClass.url_func, var0=2, var1=7)
        == "url_func?var0=2&var1=7"
    Raises error if argument names are not correct.
    """
    # Find variable names of url_func.
    default = []
    non_default = []
    sig = inspect.signature(url_func)
    for param in sig.parameters.values():
        if param.default is param.empty:
            non_default.append(param.name)
        else:
            default.append(param.name)
    # Check if variable names are correct.
    for param in args:
        if param not in default + non_default:
            raise ValueError(
                f"'{param}' is not a valid Parameter of {url_func.__name__}."
            )
    url = url_func.__name__
    if args:
        # If args are supplied, enforce that all mandatory variables are
        # contained in args.
        for param in non_default:
            if param not in args and param != "self":
                raise ValueError(
                    f"Parameter '{param}' must be supplied."
                )
        url += "?" + "&".join(
            [f"{key}={value}" for key, value in args.items()]
        )
    return url

def convert_html_to_pdf(source_html, output_file):
    """Create PDF from HTML-string."""
    with open(output_file, "w+b") as f:
        xhtml2pdf.pisa.CreatePDF(source_html, dest=f)

def date_time_string_now():
    """Return current date and time as a string to create time stamps."""
    now = datetime.datetime.now()
    return now.strftime("%Y%m%d_%H%M%S")

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
    # Sort based on modification date
    runs.sort(key=lambda x: (x[1], x[0]), reverse=True)
    # Remove date after sorting
    return [x[0] for x in runs]

def get_all_results():
    """Return list of all analysis result files in report directory sorted
    by modification time.
    """
    files = []
    for f in os.listdir(NANODIP_REPORTS):
        for e in RESULT_ENDING.values():
            if f.endswith(e):
                mod_time = os.path.getmtime(
                    os.path.join(NANODIP_REPORTS, f)
                )
                files.append([f, mod_time])
    files.sort(key=lambda x: (x[1], x[0]), reverse=True)
    return [f[0] for f in files]

def files_by_ending(directory, sample_name, ending):
    """Returns a list containing all sample output files with a given
    ending.
    """
    sample_path = os.path.join(directory, sample_name)
    output_files = []
    for root, _, files in os.walk(sample_path):
        output_files.extend(
            [os.path.join(root, f) for f in files if f.endswith(ending)]
        )
    return output_files

def predominant_barcode(sample_name):
    """Returns the predominant barcode within all fast5 files."""
    fast5_files = files_by_ending(DATA, sample_name, ending=".fast5")
    pass_fast5_files = [f for f in fast5_files if "_pass_" in f]
    barcode_hits = []
    for barcode in BARCODE_NAMES:
        barcode_hits.append(
            len([f for f in pass_fast5_files if barcode in f])
        )
    max_barcode_cnt = max(barcode_hits)
    if max_barcode_cnt > 1:
        predominant = BARCODE_NAMES[
            barcode_hits.index(max_barcode_cnt)
        ]
    else:
        predominant = "undetermined"
    return predominant

def reference_annotations():
    """Return list of all reference annotation files (MS Excel XLSX format)."""
    annotations = []
    for r in os.listdir(ANNOTATIONS):
        if r.endswith(".xlsx"):
            annotations.append(r)
    return [a.replace(".xlsx", "") for a in annotations]

def composite_path(directory, *args):
    """Generate composite file-paths of the type
        'directory/arg1_arg2_arg3'
    """
    file_name = "_".join([str(x) for x in args])
    return os.path.join(
        directory,
        file_name,
    )

def discrete_colors(names):
    """Pseudorandom color scheme based on hashed values. Colors
    of methylation classes will be fixed to their name.
        Args:
            names: List of strings.
        Returns:
            Dictionary of color scheme for all string elements.
    """
    color = {}
    for var in set(names):
        hash_str = hashlib.md5(bytes(var, "utf-8")).digest()
        hash1 = int.from_bytes(hash_str[:8], byteorder="big")
        hash2 = int.from_bytes(hash_str[8:12], byteorder="big")
        hash3 = int.from_bytes(hash_str[12:], byteorder="big")
        hue = hash1 % 365
        saturation = hash2 % 91 + 10
        lightness = hash3 % 41 + 30
        # hsl has to be transformed to rgb, since otherwise not all colors
        # are displayed correctly, probably due to plotly bug.
        rgb_frac = colorsys.hls_to_rgb(hue/364, lightness/100, saturation/100)
        rgb = tuple(int(255 * x) for x in rgb_frac)
        color[var] = f"rgb{rgb}"
    return color
