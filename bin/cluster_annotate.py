
import decoupler as dc   
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os


# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run Clustering and annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
###############################

adata = sc.read_h5ad(input_path)
print(adata.raw)

sc.pp.neighbors(adata)
sc.tl.leiden(adata)
print(adata.raw)

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata)
print(adata.raw)
sc.pl.umap(
    adata, color=["leiden", "sample_id"], palette=sc.pl.palettes.default_20
)

adata = adata.raw.to_adata()
df_markers = pd.read_csv("../data/Stroma_oldclusters_DEgenes_Venice_top100xlsx.csv", sep=";")
marker_genes = []

for ind, row in df_markers.iterrows():

    for cell_type in df_markers.columns:
        marker_genes.append([row[cell_type], cell_type])


marker_genes = pd.DataFrame(marker_genes, columns=['gene', 'cell_type']) 
marker_genes['weight'] = 1
print(marker_genes)

dc.run_ora(adata, marker_genes, min_n=0, source='cell_type', target='gene', weight='weight')

print(adata)

tmp = dc.get_acts(adata, obsm_key='ora_estimate')
print(adata.obsm)
print("============= ORA p-Vals ============= ")
print(adata.obsm["ora_pvals"])
print("============= ora_estimate ============= ")
print(adata.obsm["ora_estimate"])

sc.pl.umap(tmp, color=list(tmp.var.index)+['leiden'])

# Write to file
adata.write(os.path.join(output_path, 'integrated.h5ad'))