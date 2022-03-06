"""
## Utils

Contains utility functions.
"""


# start_external_modules
import datetime
import jinja2
import os
import pandas as pd
import xhtml2pdf.pisa
# end_external_modules

# start_internal_modules
from config import ILUMINA_CG_MAP
# end_internal_modules

def extract_referenced_cpgs(sample_methylation,
                            output_overlap,
                            output_overlap_cnt):
    """Extract ilumina CpG sites including methylation status from sample.
    Args:
        sample_methylation: methylation file of sample
        output_overlap: file path of CpG overlap
        output_overlap_cnt: file path of CpG overlap count
    """
    reference_cpgs = pd.read_csv(
        ILUMINA_CG_MAP,
        delimiter="\t",
        names=["ilmnid","chromosome","strand","start"],
    )
    sample_cpgs = pd.read_csv(
        sample_methylation,
        delimiter="\t",
    )
    cpgs = pd.merge(sample_cpgs, reference_cpgs, on=["chromosome", "start"])
    # Extract singelton CpG's
    cpgs = cpgs.loc[cpgs["num_cpgs_in_group"] == 1]
    cpgs = cpgs.loc[
       (~cpgs["chromosome"].isin(["chrX", "chrY"])) # TODO is this necessary?
       & (~cpgs["ilmnid"].duplicated())
    ]
    cpgs["is_methylated"] = 0
    cpgs.loc[cpgs["methylated_frequency"] > 0.5 ,"is_methylated"] = 1
    # Write overlap Data Frame
    cpgs[["ilmnid", "is_methylated"]].to_csv(
        output_overlap, header=False, index=False, sep="\t")
    # Write number of CpG's
    with open(output_overlap_cnt, "w") as f:
        f.write(f"{len(cpgs)}")

def render_template(template_name, **context):
    loader = jinja2.FileSystemLoader("templates")
    template = jinja2.Environment(
        loader=loader).get_template(template_name)
    return template.render(context)

def convert_html_to_pdf(source_html, output_file):
    """Create PDF from html-string."""
    with open(output_file, "w+b") as f:
        pisa_status = xhtml2pdf.pisa.CreatePDF(source_html, dest=f)
    return pisa_status.err

def date_time_string_now():
    """Return current date and time as a string to create timestamps."""
    now = datetime.datetime.now()
    return now.strftime("%Y%m%d_%H%M%S")

