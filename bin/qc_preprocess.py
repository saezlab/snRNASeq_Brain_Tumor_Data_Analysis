from genericpath import sameopenfile
import os
import utils
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy.external as sce
from anndata._core.anndata import AnnData
from utils import out_data_path, plot_path, data_path
import os
import argparse
import plotting
import warnings

warnings.simplefilter("ignore")
sc.settings.verbosity = 0

data_path = "../data/b06x-g/G703/eblanco/projects/Aniello_ITCC-P4/results/count_matrices_to_share/snRNAseq/6.1.0"
out_data_path="../data/out_data/"
Path(out_data_path).mkdir(parents=True, exist_ok=True)
plot_path="../plots/qc_preprocess"
Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path

meta = utils.get_meta_data()


def filter_cells_genes(adata, sample_id, pdx_prefix=""):
    """Perform basic filtering on cells and genes

    This function takes sample id as input and performs cell/gene filtering and remove mitochondrial genes and saves the AnnData for each sample <sample_id>_<condition>_filtered.h5ad

    Args:
        sample_id (str): the name of the folder where the sample files stored


    Returns:
        anndata: return the AnnData object

    """

    row = meta[meta["sample_id"]==sample_id]

    condition = str(row["condition"].values[0])

    # adata = utils.read_raw_sample(condition)
    
    pre_filter_shape = np.shape(adata.X)


    doublet_thr = 0.2
    mt_thr = 20
    min_g = 200
    # This was 0.99 in the first results. 
    gene_qnt = 0.995

    if pdx_prefix!="":
        condition = f"{condition}_{pdx_prefix}"
    
    # calculate qc metrics
    # adata.var["MT"] = 
    adata.var["MT"] = adata.var_names.str.contains("MT-")
    adata.var["mt"] = adata.var_names.str.contains("mt-")
    for ind in adata.var["MT"].index:
            if adata.var["MT"][ind] or adata.var["mt"][ind]:
                adata.var["mt"][ind] = True
        
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    # calculate doublet scores
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)
    sce.pp.scrublet(adata, verbose=False)
    #  “total_counts”. Sum of counts for a gene.
    #  “n_genes_by_counts”. The number of genes with at least 1 count in a cell. Calculated for all cells.
    fig, axs = plt.subplots(2, 4, figsize=(30, 10))# , figsize=(100, 20))
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=axs[0][0])
    
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False, ax=axs[0][1])
    plotting.plot_mt_vs_counts(adata, axs[0][1], mt_thr=mt_thr)
    plotting.plot_ngenes_vs_counts(adata, axs[0][2], gene_thr=gene_thr)

    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', legend_fontsize="xx-large", show=False, ax=axs[0][2])
    
    plotting.plot_doublet_scores(adata, axs[0][3], doublet_thr=doublet_thr, fontsize=11)
    sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[1][0])
    sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=100, ax=axs[1][1])
    
    sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=100, ax=axs[1][2])

    sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 1000], kde=False, bins=100, ax=axs[1][3])
    plt.axvline(50, 0, min_g, linestyle='--', color="black")
    plt.axvline(100, 0, min_g, linestyle='--', color="black")
    plt.axvline(200, 0, min_g, linestyle='--', color="black")
    plt.axvline(300, 0, min_g, linestyle='--', color="black")
    fig.savefig(os.path.join(plot_path, f"basic_stats_before_filtering_{condition}.png"), dpi=80)

    print(f"====================== {condition} ====================== ")
    print(f"AnnData shape before any filtering: {np.shape(adata.X)}")
    # print(np.shape(adata.X))
    # adata.obs.doublet_score < doublet_thr
    # number og genes at each change it to 300 
    sc.pp.filter_cells(adata, min_genes=min_g)
    print(f"AnnData shape after filtering cells with threshold {min_g}: {np.shape(adata.X)}")
    

    # drop this to three
    sc.pp.filter_genes(adata, min_cells=3)
    # print(f"AnnData after shape filter genes: {condition}")
    # print(np.shape(adata.X))
    # filter based on total counts
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)
    adata = adata[adata.obs.pct_counts_mt < mt_thr, :]
    # print(f"AnnData after shape filter cells: {condition}")
    # print(np.shape(adata.X))

    adata = adata[adata.obs.doublet_score < doublet_thr, : ]
    # print(f"AnnData after doublet: {condition}")
    # print(np.shape(adata.X))
    adata = adata[adata.obs.n_genes_by_counts < gene_thr, : ]
    # print(f"AnnData after n_genes_by_counts: {condition}")
    # print(np.shape(adata.X))
    post_filter_shape = np.shape(adata.X)
    # Assume they are coming from the same batch
    adata.obs["batch"] = 0
    adata.obs["condition"] = condition
    print(condition)
    print(f"{sample_id}:\nAnnData shape before filtering {pre_filter_shape}")
    print(f"AnnData shape after all filtering {post_filter_shape}")
    del adata.obs["predicted_doublet"]
    adata.write(os.path.join(out_data_path, f"{sample_id}_{condition}_filtered.h5ad"))
    # print(os.path.join(out_data_path, f"{sample_id}_{condition}_filtered.h5ad"))
    

    return adata

def create_filtered_adata_files():

    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata = utils.read_raw_sample(condition)
    
        print(sample_id, condition)
        if "pdx" in condition:
            
            mouse_mask_filtered = adata.var_names.str.startswith("mm10")
            human_mask_filtered = [not item for item in mouse_mask_filtered]

            adata_filtered_human = adata[:, human_mask_filtered]
            adata_filtered_mouse = adata[:, mouse_mask_filtered]
            # print(new_indices)# .index)
            adata_filtered_mouse.var.index = pd.Index(gen.split("mm10___")[1].upper() for gen in adata_filtered_mouse.var.index.values)
            # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])
            adata_filtered_mouse.var.gene_ids = pd.Index([gen.split("mm10___")[1].upper() for gen in adata_filtered_mouse.var.gene_ids.values])
        
            adata_filtered_human.var.index = pd.Index(gen.split("GRCh38_")[1].upper() for gen in adata_filtered_human.var.index.values)
            # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])
            adata_filtered_human.var.gene_ids = pd.Index([gen.split("GRCh38_")[1].upper() for gen in adata_filtered_human.var.gene_ids.values])

            filter_cells_genes(adata_filtered_human, sample_id, "human")
            filter_cells_genes(adata_filtered_mouse, sample_id, "mouse")
    
        else:
            adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)
            adata.var.gene_ids = pd.Index([gen.upper() for gen in adata.var.gene_ids.values])
            filter_cells_genes(adata, sample_id)
        
            
def analyze_pdx_samples(sample_id):

    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    for ind, sample_id in enumerate(["sample01", "sample03", "sample05"]):
    
        row = meta[meta["sample_id"]==sample_id]

        condition = str(row["condition"].values[0])

        adata_raw = utils.read_raw_sample(condition)
        mouse_mask_raw = adata_raw.var_names.str.startswith("mm10")
        human_mask_raw = [not item for item in mouse_mask_raw]

        adata_raw_mouse = adata_raw[:, mouse_mask_raw]
        adata_raw_human = adata_raw[:, human_mask_raw]
        

        adata_filtered = get_processed_sample_from_adata_file(sample_id)
        
        mouse_mask_filtered = adata_filtered.var_names.str.startswith("mm10")
        human_mask_filtered = [not item for item in mouse_mask_filtered]

        adata_filtered_mouse = adata_filtered[:, mouse_mask_filtered]
        adata_filtered_human = adata_filtered[:, human_mask_filtered]
        sc.pp.calculate_qc_metrics(adata_filtered_mouse, qc_vars=["mt"], inplace=True)
        sc.pp.calculate_qc_metrics(adata_filtered_human, qc_vars=["mt"], inplace=True)

        adata_filtered_mouse.write(os.path.join(out_data_path, 'mouse_separated_{condition}.h5ad'))
        adata_filtered_human.write(os.path.join(out_data_path, 'human_mouse_separated_{condition}.h5ad'))

        df = pd.concat([adata_filtered_mouse.obs.log1p_total_counts, adata_filtered_human.obs.log1p_total_counts], axis=1)

        # df = pd.concat([adata_filtered_mouse.obs.total_counts, adata_filtered_human.obs.total_counts], axis=1)
        
        df.columns = ["mouse_total_counts", "human_total_counts"]
        
        sns.violinplot(data=df, ax=axs[ind], title=condition)

    fig.savefig(os.path.join(plot_path, f"human_vs_mouse_violin_plot.png") , dpi=300)
    
def separate_pdx_samples():

    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    for ind, sample_id in enumerate(["sample01", "sample03", "sample05"]):
    
        row = meta[meta["sample_id"]==sample_id]

        condition = str(row["condition"].values[0])

        adata_raw = utils.read_raw_sample(condition)
        mouse_mask_raw = adata_raw.var_names.str.startswith("mm10")
        human_mask_raw = [not item for item in mouse_mask_raw]

        adata_raw_mouse = adata_raw[:, mouse_mask_raw]
        adata_raw_human = adata_raw[:, human_mask_raw]
        

        adata_filtered = get_processed_sample_from_adata_file(sample_id)
        
        mouse_mask_filtered = adata_filtered.var_names.str.startswith("mm10")
        human_mask_filtered = [not item for item in mouse_mask_filtered]

        adata_filtered_mouse = adata_filtered[:, mouse_mask_filtered]
        adata_filtered_human = adata_filtered[:, human_mask_filtered]
        sc.pp.calculate_qc_metrics(adata_filtered_mouse, qc_vars=["mt"], inplace=True)
        sc.pp.calculate_qc_metrics(adata_filtered_human, qc_vars=["mt"], inplace=True)

        adata_filtered_mouse.write(os.path.join(out_data_path, 'mouse_separated_{condition}.h5ad'))
        adata_filtered_human.write(os.path.join(out_data_path, 'human_mouse_separated_{condition}.h5ad'))

        df = pd.concat([adata_filtered_mouse.obs.log1p_total_counts, adata_filtered_human.obs.log1p_total_counts], axis=1)

        # df = pd.concat([adata_filtered_mouse.obs.total_counts, adata_filtered_human.obs.total_counts], axis=1)
        
        df.columns = ["mouse_total_counts", "human_total_counts"]
        
        sns.violinplot(data=df, ax=axs[ind], title=condition)

    fig.savefig(os.path.join(plot_path, f"human_vs_mouse_violin_plot.png") , dpi=300)



def get_processed_sample_from_adata_file(sample_id):
    """Given samples id get filtered adata object

    This function takes sample id as input and returns the filtered AnnData object

    Args:
        sample_id (str): the name of the folder where the sample files stored
    
    Returns:
        Filtered AnnData object of the sample

    """
    row = meta[meta["sample_id"]==sample_id]
    condition = str(row["condition"].values[0])
    adata = sc.read(os.path.join(out_data_path, f"{sample_id}_{condition}_filtered.h5ad"))

    return adata



if __name__ == "__main__":
    meta = utils.get_meta_data()

    create_filtered_adata_files()
