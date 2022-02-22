import argparse
import scanpy as sc
import decoupler as dc
import scanpy as sc
import decoupler as dc
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

plot_path="../plots/"
sc.settings.figdir = plot_path


# Read merged object
adata = sc.read_h5ad(input_path)

model = dc.get_progeny(organism='mouse', top=500)

model['target'] = model['target'].str.upper()


print(model)

dc.run_mlm(mat=adata, net=model, source='source', target='target', weight='weight', verbose=True)


print(adata.obsm['mlm_estimate'])

acts = dc.get_acts(adata, obsm_key='mlm_estimate')
acts

sc.pl.umap(acts, color=list(set(model["source"]))+["sample_id", "leiden"], vcenter=0, cmap='coolwarm', save=f"{sample_type}_umap_progeny")


mean_acts_cluster = dc.summarize_acts(acts, groupby='leiden', min_std=0)

ax = sns.clustermap(mean_acts_cluster, xticklabels=mean_acts_cluster.columns, vmin=-2, vmax=2, cmap='coolwarm',)
ax.savefig(f"{plot_path}/{sample_type}_clustermap_cluster_progeny", dpi=300)
plt.clf()



mean_acts_sampid = dc.summarize_acts(acts, groupby='sample_id', min_std=0)

ax = sns.clustermap(mean_acts_sampid, xticklabels=mean_acts_sampid.columns, vmin=-2, vmax=2, cmap='coolwarm',)
ax.savefig(f"{plot_path}/{sample_type}_clustermap_sample_id_progeny", dpi=300)
plt.clf()




# python downstream_analysis.py -i ../data/out_data/mouse_integrated.h5ad -o ../data/out_data -st mouse
