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

plot_path="../plots/downstream"
Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path

# Read merged object
adata = sc.read_h5ad(input_path)

organism = sample_type
model = dc.get_progeny(organism=sample_type, top=500)

model['target'] = model['target'].str.upper()

net = dc.get_dorothea(organism=sample_type, levels=['A','B','C'])
net

dc.run_mlm(mat=adata, net=net, source='source', target='target', weight='weight', verbose=True)

print(adata.obsm['mlm_estimate'])

acts = dc.get_acts(adata, obsm_key='mlm_estimate')
acts

mean_acts = dc.summarize_acts(acts, groupby='final_annotation', min_std=0.75)
mean_acts

sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-5, vmax=5, cmap='coolwarm')
plt.show()


sc.pl.umap(acts, color=['AHR', "GATA2", "GATA3", "FOXA1", "FOXP1", "SOX2", "AHR", "ZFX", 'final_annotation'], cmap='coolwarm', vcenter=0)
plt.show()
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
