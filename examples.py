import data, plots, config
from data import ReferenceData, SampleData
from plots import CNVData, UMAPData

sample_name = "B2021_48700_20211112_BC11"
reference_name = "GSE90496_IfP01"
sample = data.SampleData(sample_name)

#data.make_binary_reference_data()

ref = ReferenceData(reference_name) 
#cnv = CNVData(sample)
#umap = UMAPData.from_name(sample_name, reference_name)


