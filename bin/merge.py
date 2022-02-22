# import imp
import scanpy as sc
import scanpy.external as sce
from utils import get_meta_data
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import os
import argparse

'''
Open all samples QC processed files, merge them
'''
plot_path="../plots/"
sc.settings.figdir = plot_path
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run Merging')
parser.add_argument('-i', '--input_dir', help='Input directory containing all sample directories', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed object', required=True)
parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
args = vars(parser.parse_args())

input_path = args['input_dir']
output_path = args['output_dir']
sample_type = args['sample_type']
###############################

# Load meta data
meta = get_meta_data()
samples = np.unique(meta['sample_id'])
print(samples)


adata = []
for sample in os.listdir(input_path):
    if sample_type=="mouse" and sample_type in sample and "sample05" not in sample:
        # Read adata
        print(sample)
        tmp = sc.read_h5ad(os.path.join(input_path, sample))
        print(tmp)
        # Fetch sample metadata
        tmp.obs["sample_id"] = sample.split(".")[0]
    
        # Append
        adata.append(tmp)
        del tmp
    
    elif sample_type=="human" and (sample_type in sample or "tumor" in sample) and "sample01" not in sample:
        print(sample)
        tmp = sc.read_h5ad(os.path.join(input_path, sample))
        print(tmp)
        tmp.obs["sample_id"] = sample.split(".")[0]
        adata.append(tmp)
        del tmp
    
# Merge objects and delete list
adata = adata[0].concatenate(adata[1:], join='outer')


# Log-normalize expression
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)





# Compute HVG
sc.pp.highly_variable_genes(adata, batch_key='batch')
sc.pl.highly_variable_genes(adata, show=False, save=f'{sample_type}_merged_hvg.pdf')


adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata.raw = adata


# Filter by HVG
num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

adata.var["MT"] = adata.var_names.str.contains("MT-")
adata.var["mt"] = adata.var_names.str.contains("mt-")
for ind in adata.var["MT"].index:
        if adata.var["MT"][ind] or adata.var["mt"][ind]:
            adata.var["mt"][ind] = True

# Update QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)



# Run PCA
sc.pp.scale(adata)

# print(type(adata.X.todense()))
# ax = sns.heatmap(adata.X)
# plt.show()

sc.tl.pca(adata, svd_solver='arpack')

# Run UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(
    adata, color=["condition"], palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_merged_condition.pdf'
)

# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_merged.h5ad'))


# python merge.py -i ../data/out_data -o ../data/out_data  -st mouse
# python merge.py -i ../data/out_data -o ../data/out_data  -st human