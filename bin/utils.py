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
