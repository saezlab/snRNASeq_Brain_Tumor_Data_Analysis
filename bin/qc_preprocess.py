import scanpy as sc
import scanpy as sc
import scanpy.external as sce

import numpy as np
import pandas as pd

import os
import argparse

dataset_path = "/Users/ahmet/Google Drive/Projects/saezlab/stefan_brain_tumor_sc/b06x-g/G703/eblanco/projects/Aniello_ITCC-P4/results/count_matrices_to_share/snRNAseq/6.1.0"
sample_path = "/Users/ahmet/Google Drive/Projects/saezlab/stefan_brain_tumor_sc/b06x-g/G703/eblanco/projects/Aniello_ITCC-P4/results/count_matrices_to_share/snRNAseq/6.1.0/OE0290_pedBrainPDX_B062_007/pdx01-x"

adata = sc.read_10x_mtx(os.path.join(sample_path, 'filtered_feature_bc_matrix'), cache=True)
print(adata.X.shape)
print(adata.var["gene_ids"])
print(adata.var_names)
adata.var_names_make_unique()
print(adata.X.shape)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata.X.shape)


adata.var['mt'] = adata.var_names.str.startswith('mt-')
print(type(adata.var['mt']))
for item in adata.var['mt']:
    if item:
        
        print(item)
        