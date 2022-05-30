setwd("/Users/ahmet/Google Drive/Projects/saezlab/snRNASeq_Brain_Tumor_Data_Analysis/data/out_data")


test = readH5AD("sample01_B062_007+pdx01-x_human_filtered.h5ad", verbose = TRUE)

as.Seurat(test, counts = "X", data = NULL)



mouse_merged = readH5AD("mouse_merged.h5ad", verbose = TRUE)
mouse_merged_Seurat = as.Seurat(mouse_merged, counts = "X", data = NULL)

human_merged = readH5AD("human_merged.h5ad", verbose = TRUE)
human_merged_Seurat = as.Seurat(human_merged, counts = "X", data = NULL)


setwd("/Users/ahmet/Google Drive/Projects/saezlab/snRNASeq_Brain_Tumor_Data_Analysis/data/out_data")
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(zellkonverter)



human_integrated = readH5AD("human_integrated_clustered.h5ad", verbose = TRUE)
human_integrated_Seurat = as.Seurat(human_integrated, counts = "X", data = NULL)
SaveH5Seurat(human_integrated_Seurat, "human_integrated_Seurat.h5Seurat" , overwrite = TRUE)

human_integrated_raw = readH5AD("human_integrated_clustered_raw.h5ad", verbose = TRUE)
human_integrated_Seurat_raw = as.Seurat(human_integrated_raw, counts = "X", data = NULL)
SaveH5Seurat(human_integrated_Seurat_raw, "human_integrated_raw.h5Seurat" , overwrite = TRUE)

human_integrated_log1pnormalized = readH5AD("human_integrated_clustered_log1pnormalized.h5ad", verbose = TRUE)
human_integrated_Seurat_log1pnormalized = as.Seurat(human_integrated_log1pnormalized, counts = "X", data = NULL)
SaveH5Seurat(human_integrated_Seurat_log1pnormalized, "human_integrated_log1pnormalized.h5Seurat" , overwrite = TRUE)

mouse_integrated_clustered = readH5AD("mouse_integrated_clustered.h5ad", verbose = TRUE)
mouse_integrated_Seurat_clustered = as.Seurat(mouse_integrated_clustered, counts = "X", data = NULL)
SaveH5Seurat(mouse_integrated_Seurat_clustered, "mouse_integrated_Seurat.h5Seurat" , overwrite = TRUE)

mouse_integrated_raw = readH5AD("mouse_integrated_clustered_raw.h5ad", verbose = TRUE)
mouse_integrated_Seurat_raw = as.Seurat(mouse_integrated_raw, counts = "X", data = NULL)
SaveH5Seurat(mouse_integrated_Seurat_raw, "mouse_integrated_raw.h5Seurat" , overwrite = TRUE)

mouse_integrated_log1pnormalized = readH5AD("mouse_integrated_clustered_log1pnormalized.h5ad", verbose = TRUE)
mouse_integrated_Seurat_log1pnormalized = as.Seurat(mouse_integrated_log1pnormalized, counts = "X", data = NULL)
SaveH5Seurat(mouse_integrated_Seurat_log1pnormalized, "mouse_integrated_log1pnormalized.h5Seurat" , overwrite = TRUE)

human_integrated = readH5AD("human_integrated.h5ad", verbose = TRUE)
human_integrated_Seurat = as.Seurat(human_integrated, counts = "X", data = NULL)
SaveH5Seurat(human_integrated_Seurat, "human_integrated_Seurat.h5Seurat" , overwrite = TRUE)



integrate_test = readh5AD("mouse_integrated.h5ad", verbose=TRUE)

integrate_Seurat = as.Seurat(integrate_test, counts = "X", data = NULL)


sample_name_list = c("sample01_B062_007+pdx01-x_human_filtered", "sample08_mouse_NSG+striatum01_filtered",
"sample07_mouse_NSG+cerebellum01_filtered",
"sample06_B062_025+tumor01_filtered",
"sample05_B062_025+pdx01-x_mouse_filtered",
"sample05_B062_025+pdx01-x_human_filtered",
"sample04_B062_013+tumor01_filtered",
"sample03_B062_013+pdx01-x_mouse_filtered",
"sample03_B062_013+pdx01-x_human_filtered",
"sample02_B062_007+tumor01_filtered",
"sample01_B062_007+pdx01-x_mouse_filtered")

for (val in sample_name_list)
{
    print(paste0(val, ".h5ad"))
    human_integrated = readH5AD(paste0(val, ".h5ad"), verbose = TRUE)
    human_integrated_Seurat = as.Seurat(human_integrated, counts = "X", data = NULL)
    SaveH5Seurat(human_integrated_Seurat, paste0(val, ".h5Seurat"), overwrite = TRUE)
    
}