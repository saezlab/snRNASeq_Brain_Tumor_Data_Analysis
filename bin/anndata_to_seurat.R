

int_data <- load("Riemondy_Seurat_obkject_TME.RData")
write.csv(sce.m@reductions$umap@cell.embeddings, "cell_embeddings.csv", quote = FALSE)
write.csv(sce.m@meta.data,"meta_data.csv", quote = FALSE)
write.csv(sce.m@assays$RNA@counts,"count.csv", quote = FALSE)


setwd("../data/out_data")


test = readH5AD("sample01_B062_007+pdx01-x_human_filtered.h5ad", verbose = TRUE)

as.Seurat(test, counts = "X", data = NULL)



mouse_merged = readH5AD("mouse_merged.h5ad", verbose = TRUE)
mouse_merged_Seurat = as.Seurat(mouse_merged, counts = "X", data = NULL)

human_merged = readH5AD("human_merged.h5ad", verbose = TRUE)
human_merged_Seurat = as.Seurat(human_merged, counts = "X", data = NULL)


setwd("../data/out_data")
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(zellkonverter)



human_integrated = readH5AD("human_integrated_clustered.h5ad", verbose = TRUE)
human_integrated_Seurat = as.Seurat(human_integrated, counts = "X", data = NULL)
SaveH5Seurat(human_integrated_Seurat, "human_integrated_Seurat.h5Seurat" , overwrite = TRUE)

---
title: "test"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r}

    library("reticulate")
    library("ggplot2")
    library("SingleCellExperiment")
    library("scater")
    library("Seurat")
```
```{r}
setwd("../data/aniello_processed_objects/")
gojo_ss2_data <- load("Gojo_SS2_Seurat.RData")
gojo_ss2_data <- sce.m
gojo_meta <- load("Gojo_SS2_Seurat_metadata.RData")
gojo_meta <- meta
gojo_ss2_data$final_annotation <- gojo_meta$final_annotation
```

```{r}
write.table(gojo_ss2_data$final_annotation, "cell_types.csv", sep = ',', row.names = T, col.names = T, quote = F)
```


```{python}
import scanpy as sc
import pandas as pd
import os 

os.chdir("../data/aniello_processed_objects/")
df = pd.read_csv("counts.csv")
```

```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = '../data/aniello_processed_objects/')
```

```{python}
import os 

os.getcwd()
df_umap = pd.read_csv("umap.csv")
df_cell_types = pd.read_csv("cell_types.csv")

adata = sc.AnnData(df.iloc[:,:].T, obs = pd.DataFrame(index=list(df.columns)), var = pd.DataFrame(index=list(df.index)))
adata.obsm["X_umap"] = df_umap
adata.obs["final_annotation"] = df_cell_types
adata.write("Gojo_SS2.h5ad")
```




## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
