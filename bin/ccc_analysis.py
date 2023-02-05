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
import liana as li
import plotnine as p9
# Only needed for visualization:
import matplotlib.pyplot as plt
# import liana's rank_aggregate
from plotnine.scales import scale_x_continuous, scale_x_discrete
from liana.mt import rank_aggregate
from plotnine import ggplot, geom_point, aes, theme, element_text, facet_grid



# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run CCC Analysis')
parser.add_argument('-i', '--input_path', help='Input path to processed object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
sample_type = args['sample_type']






sample_name = "prPDX_human"

if input_path.endswith("Riemondy_TME.h5ad"):
    sample_name = "Riemondy_TME"
elif input_path.endswith("Gojo_SS2_updated.h5ad"):
    sample_name = "Gojo_SS2"

print("Sample name:", sample_name)
plot_path = f"../plots/downstream/CCC/{sample_name}"

Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path

"""# Read merged object
adata = sc.read_h5ad(input_path)
if not adata.raw:
    print("Assume raw matrix and perform normalization...")
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)"""


# print(adata.obs)
organism = sample_type

print(li.mt.show_methods())

print("resources", li.resource.show_resources())


levels = ['TME_level_1', 'TME_level_2', 'TME_level_3']



if sample_name=="prPDX_human":
    levels = ["TME_level_1", "TME_level_2"]
elif sample_name=="Riemondy_TME":
    levels = ["TME_level_1", "TME_level_2", "TME_level_3"]
levels = ["TME_level_3"]
p9.options.dpi = 600
for level in levels:
    print(level)
    adata = sc.read_h5ad(f"/Users/ahmet/Downloads/{sample_name}_{level}.h5ad")

    malignant_cat = []
    normal_cat = []

    for cat in adata.obs[level].cat.categories:
        if cat.startswith("Malignant"):
            malignant_cat.append(cat)
        else:
            normal_cat.append(cat)


    """if not adata.raw:  
        li.mt.rank_aggregate(adata, use_raw=False, groupby=level, expr_prop=0.1, verbose=True)
    else:
        li.mt.rank_aggregate(adata, groupby=level, expr_prop=0.1, verbose=True)

    adata.write(f"/Users/ahmet/Downloads/{sample_name}_{level}.h5ad")"""
    my_p = li.pl.dotplot(adata = adata,
                        colour='lrscore',
                        size='spec_weight',
                        # size='pvals',
                        inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                        source_labels=malignant_cat,
                        target_labels=normal_cat,
                        top_n=25,
                        orderby='steady_rank',
                        orderby_ascending=True,
                        # finally, since cpdbv2 suggests using a filter to FPs
                        # we filter the pvals column to <= 0.05
                        filterby='cellphone_pvals',
                        filter_lambda=lambda x: x <= 0.05,
                        figure_size = (144, 8)
                        )


    
    my_p = (my_p +
    facet_grid('~source', scales="free_x", space='free_x') +
    p9.theme(strip_text = p9.element_text(size = 6, face="bold", colour = "gray"), 
        axis_text_x=element_text(size=7, face="bold", angle=90))
    )

    my_p.save(f"/Users/ahmet/Desktop/source_malignant_target_normal_{sample_name}_{level}.pdf", limitsize = False)

    
    my_p = li.pl.dotplot(adata = adata,
                        colour='lrscore',
                        size='spec_weight',
                        # size='pvals',
                        inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                        source_labels=normal_cat,
                        target_labels=malignant_cat,
                        top_n=25,
                        orderby='steady_rank',
                        orderby_ascending=True,
                        # finally, since cpdbv2 suggests using a filter to FPs
                        # we filter the pvals column to <= 0.05
                        filterby='cellphone_pvals',
                        filter_lambda=lambda x: x <= 0.05,
                        figure_size = (96, 8)
                        )


    my_p = (my_p +
    facet_grid('~source', scales="free_x", space='free_x') +
    p9.theme(strip_text = p9.element_text(size = 8, face="bold", colour = "gray"), 
        axis_text_x=element_text(size=7, face="bold", angle=90))
    )
    


    my_p.save(f"/Users/ahmet/Desktop/source_normal_target_malignant_{sample_name}_{level}.pdf", limitsize = False)



# python ccc_analysis.py -i ../data/aniello_processed_objects/sce_updated.h5ad -o . -st human
# python ccc_analysis.py -i ../data/aniello_processed_objects/Gojo_SS2_updated.h5ad -o . -st human
# python ccc_analysis.py -i ../data/aniello_processed_objects/Riemondy_TME.h5ad -o . -st human
# python ccc_analysis.py -i ../data/aniello_processed_objects/sce_updated.h5ad -o . -st human
# python ccc_analysis.py -i /Users/ahmet/Downloads/Gojo_SS2_updated.h5ad -o . -st human
# python ccc_analysis.py -i /Users/ahmet/Downloads/sce_updated.h5ad -o . -st human
# python ccc_analysis.py -i ../data/aniello_processed_objects/Riemondy_TME.h5ad -o . -st human
# python ccc_analysis.py -i /Users/ahmet/Downloads/Riemondy_TME.h5ad -o . -st human

"""python ccc_analysis.py -i /Users/ahmet/Downloads/Riemondy_TME.h5ad -o . -st human
python ccc_analysis.py -i /Users/ahmet/Downloads/sce_updated.h5ad -o . -st human
python ccc_analysis.py -i /Users/ahmet/Downloads/Gojo_SS2_updated.h5ad -o . -st human"""

