import data
import os
import pandas as pd
import config
import jinja2

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
        config.ILUMINA_CG_MAP,
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

