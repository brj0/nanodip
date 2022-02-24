from data import ReferenceData, SampleData, ReferenceGenome
from plots import CNVData, UMAPData
from tqdm import tqdm
import config
import csv
import data, plots, config
import numpy as np
import openpyxl
import os
import pandas as pd
import pysam
import re
import time

sample_name = "B2021_48700_20211112_BC11"
reference_name = "GSE90496_IfP01"
reference_name = "20210721_EpiDiP_anno"
sample = data.SampleData(sample_name)

#data.make_binary_reference_data()

#ref = ReferenceData(reference_name) 
cnv = CNVData(sample_name)
#umap = UMAPData.from_name(sample_name, reference_name)

g = ReferenceGenome()


