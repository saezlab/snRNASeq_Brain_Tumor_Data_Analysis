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
# parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
# sample_type = args['sample_type']

plot_path="../plots/downstream"
Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path

# Read merged object
adata = sc.read_h5ad(input_path)

print(adata)

sc.pl.umap(adata, color="final_annotation")