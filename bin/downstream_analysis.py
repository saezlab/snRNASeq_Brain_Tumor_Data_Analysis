

import os
import argparse
import scanpy as sc
import decoupler as dc
import scanpy as sc
import decoupler as dc
from pathlib import Path
import scanpy.external as sce
# Only needed for visualization:
import matplotlib.pyplot as plt
import seaborn as sns
import sys



# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run Downstream Analysis')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
sample_type = args['sample_type']

sample_name = "prPDX_human"
sample_name = "Riemondy_TME"
plot_path = f"../plots/downstream/TF/{sample_name}"

Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path

# Read merged object
adata = sc.read_h5ad(input_path)
organism = sample_type


net = dc.get_dorothea(organism=sample_type, levels=['A','B','C'])

net['target'] = net['target'].str.upper()

# dc.run_mlm(mat=adata, use_raw=False, net=net, source='source', target='target', weight='weight', verbose=True)
if adata.raw:
    dc.run_mlm(mat=adata, net=net, source='source', target='target', weight='weight', verbose=True)
else:
    # Assume raw matrix and perform normalization
    print("Assume raw matrix and perform normalization...")
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    dc.run_mlm(mat=adata, use_raw=False, net=net, source='source', target='target', weight='weight', verbose=True)

#  print(adata.obsm['mlm_estimate'])

acts = dc.get_acts(adata, obsm_key='mlm_estimate')
# print("Acts:")
# print(acts)


meta_columns = ["TME_level_1", "TME_level_2", "TME_level_3", "final_annotation"]

if sample_name=="prPDX_human":
    meta_columns = ["TME_level_1", "TME_level_2", "final_annotation"]
elif sample_name=="Riemondy_TME":
    meta_columns = ["TME_level_1", "TME_level_2", "TME_level_3"]

for col in meta_columns:
    if col in adata.obs.columns:
        mean_acts = dc.summarize_acts(acts, groupby=col, min_std=0.75)
        # print("Mean acts:", mean_acts)
        sns.clustermap(mean_acts, xticklabels=mean_acts.columns, figsize=(20, 10), vmin=-5, vmax=5, cmap='coolwarm', annot_kws={"size": 4})
        plt.savefig(f"{plot_path}/{sample_type}_{col}_dorothea_mean_activities.pdf")


sc.pl.umap(acts, color=list(mean_acts.columns)+meta_columns, show=False, cmap='coolwarm', vcenter=0, save=f"{sample_type}_umap_dorothea_all_TFs")

for tf in list(mean_acts.columns)+meta_columns:
    sc.pl.umap(acts, color=tf, show=False, cmap='coolwarm', vcenter=0, save=f"{sample_type}_umap_dorothea_{tf}")

"""leiden_res_params = [0.1, 0.2, 0.5, 0.7, 1.0]



organism = sample_type
model = dc.get_progeny(organism=sample_type, top=100)

model['target'] = model['target'].str.upper()

print(model)

dc.run_mlm(mat=adata, net=model, source='source', target='target', weight='weight', verbose=True)

print(adata.obsm['mlm_estimate'])"""

"""
pathway_list = list(set(model["source"]))

for path in pathway_list:
    adata.obs[path] = adata.obsm['mlm_estimate'][path]
"""
"""acts = dc.get_acts(adata, obsm_key='mlm_estimate')
print(acts)


for l_param in leiden_res_params:

    sc.pl.umap(acts, color=list(set(model["source"]))+[f"leiden_{l_param}"], vcenter=0, show=False, cmap='coolwarm', save=f"{sample_type}_umap_progeny_res_{l_param}")
    mean_acts_cluster = dc.summarize_acts(acts, groupby=f"leiden_{l_param}", min_std=0)

    ax = sns.clustermap(mean_acts_cluster, xticklabels=mean_acts_cluster.columns, vmin=-2, vmax=2, cmap='coolwarm',)
    ax.savefig(f"{plot_path}/{sample_type}_clustermap_cluster_progeny_res_{l_param}.pdf", show=False, dpi=300)
    plt.clf()

    mean_acts_sampid = dc.summarize_acts(acts, groupby='sample_id', min_std=0)

    ax = sns.clustermap(mean_acts_sampid, xticklabels=mean_acts_sampid.columns, vmin=-2, vmax=2, cmap='coolwarm',)
    ax.savefig(f"{plot_path}/{sample_type}_clustermap_sample_id_progeny_res_{l_param}.pdf", show=False, dpi=300)
    plt.clf()



acts.write(os.path.join(output_path, f'{sample_type}_integrated_progeny_act.h5ad'))
acts.write(os.path.join(output_path, f'{sample_type}_integrated_progeny_act.h5ad'))"""

# python downstream_analysis.py -i ../data/out_data/mouse_integrated.h5ad -o ../data/out_data -st mouse
# python downstream_analysis.py -i ../data/out_data/human_integrated.h5ad -o ../data/out_data -st human
# python downstream_analysis.py -i ../data/aniello_processed_objects/sce_updated.h5ad -o . -st human
# python downstream_analysis.py -i ../data/aniello_processed_objects/Gojo_SS2_updated.h5ad -o . -st human
# python downstream_analysis.py -i ../data/aniello_processed_objects/Riemondy_TME.h5ad -o . -st human
