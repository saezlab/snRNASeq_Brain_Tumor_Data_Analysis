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
parser.add_argument('-res', '--resolution',  help='Resolution param for leiden', type=float, default=1.0)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
sample_type = args['sample_type']

###############################



plot_path="../plots/"
sc.settings.figdir = plot_path

adata = sc.read_h5ad(input_path)


sc.pp.neighbors(adata)

sc.tl.leiden(adata, key_added = "leiden_1.0") # default resolution in 1.0
sc.tl.leiden(adata, resolution = 0.6, key_added = "leiden_0.6")
# sc.tl.leiden(adata, resolution = 0.4, key_added = "leiden_0.4")





sc.tl.rank_genes_groups(adata, 'leiden_1.0', method='wilcoxon', key_added = "wilcoxon_1.0")
sc.tl.rank_genes_groups(adata, 'leiden_0.6', method='wilcoxon', key_added = "wilcoxon_0.6")
# sc.tl.rank_genes_groups(adata, 'leiden_0.4', method='wilcoxon', key_added = "wilcoxon_0.4")
# sc.pl.rank_genes_groups(adata)

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", show=False, groupby="leiden", show_gene_labels=True, save=f'{sample_type}_deg_clusters_heatmap')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon_1.0", show=False, groupby="leiden_1.0", save=f'{sample_type}_deg_clusters_dotplot_1.0')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon_0.6", show=False, groupby="leiden_0.6", save=f'{sample_type}_deg_clusters_dotplot_0.6')
# sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon_0.4", show=False, groupby="leiden_0.4", save=f'{sample_type}_deg_clusters_dotplot_0.4')


# sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, key="wilcoxon", show=False, groupby="leiden", save=f'{sample_type}_ranked_genes_stackedviolin')


sc.pl.umap(adata, color=['leiden_0.6', 'leiden_1.0'], legend_loc="on data", legend_fontsize="xx-small", palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_clusters_diff_res')
"""sc.pl.umap(
    adata, color=["leiden"], palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_clusters'
)"""

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

mean_enr10 = dc.summarize_acts(acts, groupby='leiden_1.0')
mean_enr06 = dc.summarize_acts(acts, groupby='leiden_0.6')
# mean_enr04 = dc.summarize_acts(acts, groupby='leiden_0.4')
# print(mean_enr)
sns.clustermap(mean_enr10, xticklabels=mean_enr10.columns)
plt.savefig(f"{plot_path}/{sample_type}_clustermap_res10.pdf" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}_res10.pdf") 
sns.clustermap(mean_enr06, xticklabels=mean_enr06.columns)
plt.savefig(f"{plot_path}/{sample_type}_clustermap_res06.pdf" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}_res06.pdf") 
"""sns.clustermap(mean_enr04, xticklabels=mean_enr04.columns)
plt.savefig(f"{plot_path}/{sample_type}_clustermap_res04.pdf" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}_res04.pdf") """

annotation_dict10 = dc.assign_groups(mean_enr10)
annotation_dict06 = dc.assign_groups(mean_enr06)
# annotation_dict04 = dc.assign_groups(mean_enr04)

# Add cell type column based on annotation
adata.obs['cell_type_1.0'] = [annotation_dict10[clust] for clust in adata.obs['leiden_1.0']]
# Add cell type column based on annotation
adata.obs['cell_type_0.6'] = [annotation_dict06[clust] for clust in adata.obs['leiden_0.6']]
# Add cell type column based on annotation
# adata.obs['cell_type_0.4'] = [annotation_dict04[clust] for clust in adata.obs['leiden_0.4']]

# Visualize
sc.pl.umap(adata, color=['cell_type_1.0', 'cell_type_0.6'], show=False, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_resxx.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')

# sc.pl.umap(tmp, color=list(tmp.var.index)+['leiden']+['sample_id'],  show=False, save=f"{sample_type}_cell_type")

# Write to file
# adata.write(os.path.join(output_path, f'{sample_type}_integrated.h5ad'))


# python cluster_annotate.py -i ../data/out_data/mouse_integrated.h5ad -o ../data/out_data -st mouse
# python cluster_annotate.py -i ../data/out_data/human_integrated.h5ad -o ../data/out_data -st human
