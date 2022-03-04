import bisect
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.io import write_json, from_json
import pysam
import random
import re
import time

from config import (
    CNV_GRID,
    CNV_URL_PREFIX,
    CNV_URL_SUFFIX,
    ENDINGS,
    NANODIP_REPORTS,
    PLOTLY_RENDER_MODE,
    UMAP_PLOT_TOP_MATCHES,
)
from data import (
    get_sample_methylation,
    SampleData,
    ReferenceData,
    ReferenceGenome,
    get_reference_methylation,
)

from plots import CNVData, UMAPData
import config, data, plots

sample_name = "B2021_48700_20211112_BC11"
reference_name = "20210721_EpiDiP_anno"
reference_name = "GSE90496_IfP01"
sample = SampleData(sample_name)
reference = ReferenceData(reference_name)
genome = ReferenceGenome()

# data.make_binary_reference_data()

# ref = ReferenceData(reference_name)
# cnv = CNVData(sample_name)
umapp = UMAPData(sample_name, reference_name)
umapp.make_umap_plot()


def umap_plot_from_data(sample, reference, umap_data_frame, close_up):
    """Create and return umap plot from UMAP data.

    Args:
        sample: sample data
        reference: reference data
        umap_data_frame: pandas data frame containing umap info. First
            row corresponds to sample.
        close_up: bool to indicate if only top matches should be plotted.
    """
    umap_sample = umap_data_frame.iloc[0]
    umap_title = f"UMAP for {sample.name} against {reference.name}, "\
        + f"{len(reference.annotated_specimens)} reference cases, "\
        + f"{len(sample.cpg_overlap)} CpGs"
    if close_up:
        umap_title = "Close-up " + umap_title
    umap_plot = px.scatter(
        umap_data_frame,
        x="x",
        y="y",
        labels={"x":"UMAP 0", "y":"UMAP 1", "color":"WHO class"},
        title=umap_title,
        color="methylation_class",
        hover_name="id",
        hover_data=["description"],
        render_mode=PLOTLY_RENDER_MODE,
        template="simple_white",
    )
    umap_plot.add_annotation(
        x=umap_sample["x"],
        y=umap_sample["y"],
        text=sample.name,
        showarrow=True,
        arrowhead=1,
    )
    umap_plot.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
        mirror=True,
    )
    umap_plot.update_xaxes(
        mirror=True,
    )
    # If close-up add hyperlinks for all references and draw circle
    if close_up:
        umap_plot.update_traces(marker=dict(size=5))
        # Add hyperlinks
        for _, row in umap_data_frame.iloc[1:].iterrows():
            umap_plot.add_annotation(
                x=row["x"],
                y=row["y"],
                text="<a href='" + CNV_URL_PREFIX + row["id"]
                    + CNV_URL_SUFFIX
                    + "' target='_blank'>&nbsp;</a>",
                showarrow=False,
                arrowhead=1,
            )
        # Draw circle
        radius = umap_data_frame["distance"].iloc[-1]
        umap_plot.add_shape(
            type="circle",
            x0=umap_sample["x"] - radius,
            y0=umap_sample["y"] - radius,
            x1=umap_sample["x"] + radius,
            y1=umap_sample["y"] + radius,
            line_color="black",
            line_width=0.5,
        )
    return umap_plot


plt=umap_plot_from_data(umapp.sample, umapp.reference, umapp.umap_df, close_up=False)
plt.show()

