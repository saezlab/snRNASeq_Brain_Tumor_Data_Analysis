from genericpath import sameopenfile
import matplotlib.pyplot as plt
import seaborn as sns 
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
parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
sample_type = args['sample_type']
###############################

plot_path="../plots/"
sc.settings.figdir = plot_path

adata = sc.read_h5ad(input_path)


sc.pp.neighbors(adata)
sc.tl.leiden(adata)


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
# sc.pl.rank_genes_groups(adata)

sc.pl.umap(
    adata, color=["leiden"], palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_clusters'
)

sc.pl.umap(
    adata, color=["sample_id"], palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_samples'
)

marker_genes = None
own_markers = True
marker_db="PanglaoDB"

if own_markers:
    # adata = adata.raw.to_adata()
    df_markers = pd.read_csv("../data/Stroma_oldclusters_DEgenes_Venice_top100xlsx.csv", sep=";")
    marker_genes = []

    # adata_filtered_mouse.var.index = pd.Index(gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.index.values)
    # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])

    for ind, row in df_markers.iterrows():

        for cell_type in df_markers.columns:
            if [row[cell_type], cell_type] not in marker_genes:
                marker_genes.append([row[cell_type], cell_type])

    marker_genes = pd.DataFrame(marker_genes, columns=['genesymbol', 'cell_type']) 
    marker_genes['weight'] = 1
else:
    
    # print(marker_genes)
    marker_genes = dc.get_resource('PanglaoDB')
    marker_db="PanglaoDB"
    # Filter by canonical_marker and human
    marker_genes = marker_genes[(marker_genes['human']=='True')&(marker_genes['canonical_marker']=='True')]

    # Remove duplicated entries
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

print(marker_genes)

marker_genes['genesymbol'] = marker_genes['genesymbol'].str.upper()

dc.run_ora(mat=adata, net=marker_genes, source='cell_type', target='genesymbol', min_n=3)

# print(adata)
# dc.run_ora(mat=adata, net=markers, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

print(adata.obsm)
print("============= ORA p-Vals ============= ")
print(adata.obsm["ora_pvals"])
print("============= ora_estimate ============= ")
print(adata.obsm["ora_estimate"])

mean_enr = dc.summarize_acts(acts, groupby='leiden')
print(mean_enr)
sns.clustermap(mean_enr, xticklabels=mean_enr.columns)

#fig = swarm_plot.get_figure()
plt.savefig(f"{plot_path}/{sample_type}_clustermap" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}") 


annotation_dict = dc.assign_groups(mean_enr)



# Add cell type column based on annotation
adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden']]

# Visualize
sc.pl.umap(adata, color='cell_type', show=False, save=f'{sample_type}_cell_type' if own_markers else f'{sample_type}_cell_type_PanglaoDB.pdf')

# sc.pl.umap(tmp, color=list(tmp.var.index)+['leiden']+['sample_id'],  show=False, save=f"{sample_type}_cell_type")

# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated.h5ad'))


# python cluster_annotate.py -i ../data/out_data/mouse_integrated.h5ad -o ../data/out_data -st mouse
# python cluster_annotate.py -i ../data/out_data/human_integrated.h5ad -o ../data/out_data -st human
