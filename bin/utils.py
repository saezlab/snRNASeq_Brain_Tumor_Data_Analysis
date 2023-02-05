

import os
import pickle
import pandas as pd
import scanpy as sc

data_path = "../data/b06x-g/G703/eblanco/projects/Aniello_ITCC-P4/results/count_matrices_to_share/snRNAseq/6.1.0"
out_data_path="../data/out_data/"
plot_path="../plots/"
sc.settings.verbosity = 0

def read_pickle(fl_path):
    """Read a pickle file

    This function reads a pickle file and returns the object.

    Args:
        fl_path (str): pickle file to be read


    Returns:
        object: return the object stored in the pickle file.

    """

    p_file = open(fl_path,'rb')
    obj_r = pickle.load(p_file)
    p_file.close()

    return obj_r

def write_pickle(fl_path, obj_w):
    """Create a pickle file given object

    This function creates a pickle file to store the input object in the given path.

    Args:
        fl_path (str): pickle file to be created
        obj_w (object): object to be stored

    """
    p_file = open(fl_path, 'wb')    
    # dump information to that file
    pickle.dump(obj_w, p_file)
    p_file.close()


def get_meta_data():
    """
    Read  meta data file and return the dataframe

    Returns:
        Meta data as data frame 
    """

    meta = pd.read_csv(f'{data_path}/meta_data.csv')
    return meta


def read_raw_sample(sample_name):
    adata = sc.read_10x_mtx(os.path.join(data_path, sample_name, 'filtered_feature_bc_matrix'), cache=True)
    adata.var_names_make_unique()
    return adata


def create_anndata_from_count_matrix(count_data_path, meta_data_path=None, cell_embedding_path=None):
    df_meta_data = pd.read_csv(meta_data_path, index_col=0, sep=";")
    df_cell_embedding = pd.read_csv(cell_embedding_path, index_col=0)
    

    # print("INDEX", df_cell_embedding.index)
    #print(df_meta_data.columns)
    #for col in df_meta_data.columns:
    #    print(col)

    
    df_count = pd.read_csv(count_data_path, index_col=0).transpose() 
    # print(df_count)
    # print(df_count.columns)
    
    adata = sc.AnnData(df_count)
    # print(adata.var_names)
    # print(adata.obs_names)
    adata.obsm['X_umap'] = df_cell_embedding.to_numpy()
    
    for col in df_meta_data.columns:
        adata.obs[col] = pd.Series(df_meta_data[col], dtype="string").values
    adata.write(f"../data/aniello_processed_objects/Riemondy_TME.h5ad")


# create_anndata_from_count_matrix("../data/aniello_processed_objects/Riemondy_Seurat_object_TME/count.csv", "../data/aniello_processed_objects/Riemondy_Seurat_object_TME/meta_data.csv", "../data/aniello_processed_objects/Riemondy_Seurat_object_TME/cell_embeddings.csv")

def add_metadata_anndata(adata_path, meta_data_path):
    adata = sc.read_h5ad(adata_path)
    print(adata.obs_names)
    df_meta_data = pd.read_csv(meta_data_path)
    print(df_meta_data)
    for col in df_meta_data.columns:
        adata.obs[col] = pd.Series(df_meta_data[col], dtype="string").values

    # adata.write("../data/aniello_processed_objects/Gojo_SS2_updated.h5ad")
    adata.write("../data/aniello_processed_objects/sce_updated.h5ad")


# add_metadata_anndata("../data/aniello_processed_objects/Gojo_SS2.h5ad", "../data/aniello_processed_objects/GojoSS2_Seurat_metadata_TME.csv")
add_metadata_anndata("../data/aniello_processed_objects/sce.h5ad", "../data/aniello_processed_objects/prPDX_human_TME_metadata.csv")
