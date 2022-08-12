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


int_data <- load("human_Seurat_final_annotation.RData")
SaveH5Seurat(sce.m, filename = "sce.h5Seurat")
Convert("sce.h5Seurat", dest = "h5ad")




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




"""malignant_B062_007+pdx01-x_human
malignant_B062_007+tumor01
malignant_B062_013+pdx01-x_human
malignant_B062_013+tumor01
malignant_B062_025+pdx01-x_human
malignant_B062_025+tumor01"""


library(tidyverse)
library(magrittr)
library(liana)
show_resources()
show_methods()
setwd("Google Drive/Projects/saezlab/snRNASeq_Brain_Tumor_Data_Analysis/data/aniello_processed_objects/")
testdata <- load("human_Seurat_final_annotation.RData")
testdata %>% dplyr::glimpse()
testdata <- sce.m
liana_test <- liana_wrap(testdata, idents_col = "final_annotation")

liana_test <- liana_test %>%
liana_aggregate()

dplyr::glimpse(liana_test)


malignant_groups <- c("malignant_B062_007+pdx01-x_human", 
"malignant_B062_007+tumor01",
"malignant_B062_013+pdx01-x_human",
"malignant_B062_013+tumor01",
"malignant_B062_025+pdx01-x_human",
"malignant_B062_025+tumor01")

c_types  = c("myeloid", "oligodendrocytes", "stromal", "lymphoid")

print(val)
liana_test %>%
dplyr::filter(source %in% c_types) %>%
dplyr::filter(target %in% malignant_groups) %>%
mutate(aggregate_rank = p.adjust(aggregate_rank, method="BH")) %>%
filter(aggregate_rank <= 0.02) %>%
liana_dotplot(source_groups = c_types,
target_groups = malignant_groups)+
theme(strip.text = element_text(size = 10, face="bold", colour = "gray6"), 
axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust=0.5), 
axis.text.y = element_text(size = 10, vjust = 0.5, hjust=0.5))+ 
scale_x_discrete(breaks=malignant_groups,
labels=c("007+pdx01-x", "007+tumor01","013+pdx01-x","013+tumor01","025+pdx01-x","025+tumor01"))
ggsave(paste0("../../plots/ccc/integrated/source-all_targets-malignant_group.pdf"))



for (val in c_types)
{
    print(val)
    liana_test %>%
    dplyr::filter(source ==val) %>%
    dplyr::filter(target %in% malignant_groups) %>%
    dplyr::top_n(25, desc(aggregate_rank)) %>%
    liana_dotplot(source_groups = c(val),
    target_groups = malignant_groups)+
    theme(strip.text = element_text(size = 10, face="bold", colour = "gray6"), 
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust=0.5), 
    axis.text.y = element_text(size = 10, vjust = 0.5, hjust=0.5))+ 
    scale_x_discrete(breaks=malignant_groups,
    labels=c("007+pdx01-x", "007+tumor01","013+pdx01-x","013+tumor01","025+pdx01-x","025+tumor01"))
    ggsave(paste0("../../plots/ccc/integrated/source-", val, "_targets-malignant_group.pdf"))
}


for (val in malignant_groups)
{
    print(val)
    liana_test %>%
    dplyr::filter(source ==val) %>%
    dplyr::top_n(25, desc(aggregate_rank)) %>%
    liana_dotplot(source_groups = c(val),
    target_groups = c_types)+
    theme(strip.text = element_text(size = 10, colour = "gray6"),
    plot.title = element_text(size=4), 
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust=0.5), 
    axis.text.y = element_text(size = 10, vjust = 0.5, hjust=0.5))
    ggsave(paste0("../../plots/ccc/integrated/source-", val, "_targets-mye_oli_stro_lymp.pdf"))
}




liana_test %>%
dplyr::filter(source =="malignant_B062_007+pdx01-x_human") %>%
dplyr::top_n(25, desc(aggregate_rank)) %>%
liana_dotplot(source_groups = c("malignant_B062_007+pdx01-x_human"),
target_groups = c("stromal", "myeloid", "oligodendrocytes"))

liana_test %>%
dplyr::filter(source =="malignant_B062_007+pdx01-x_human") %>%
dplyr::top_n(100, desc(aggregate_rank)) %>%
liana_dotplot(source_groups = c("malignant_B062_007+pdx01-x_human"),
target_groups = c("stromal", "myeloid", "oligodendrocytes"))

liana_trunc <- liana_test %>%
# only keep interactions concordant between methods
filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
heat_freq(liana_trunc)
# history()
# savehistory("~/Desktop/Untitled.Rhistory")



gojo_ss2_data <- load("Gojo_SS2_Seurat.RData")
gojo_ss2_data <- sce.m
gojo_meta <- load("Gojo_SS2_Seurat_metadata.RData")
gojo_meta <- meta
gojo_ss2_data$final_annotation <- gojo_meta$final_annotation
gojo_ss2_data %>% dplyr::glimpse()
gojo_ss2_liana <- liana_wrap(gojo_ss2_  data, idents_col = "final_annotation")
gojo_ss2_data[, gojo_ss2_data@meta.data$final_annotation %in% c("Myeloid", "")]
gojo_ss2_data@meta.data$final_annotation%>%as.factor()
gojo_ss2_liana <- gojo_ss2_liana %>%
liana_aggregate()

dplyr::glimpse(gojo_ss2_liana)


"""
“malignant_MUV006” and “myeloid”; 
“malignant_MUV006” and “lymphoid; 
“malignant_MUV006” and “stroma”.
"""


malignant_groups <- c("malignant_MUV006")

c_types  = c("myeloid", "stromal", "lymphoid")

for (val in c_types)
{
    print(val)
    gojo_ss2_liana %>%
    dplyr::filter(source ==val) %>%
    dplyr::top_n(25, desc(aggregate_rank)) %>%
    liana_dotplot(source_groups = c(val),
    target_groups = malignant_groups)+
    theme(strip.text = element_text(size = 10, face="bold", colour = "gray6"), 
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust=0.5), 
    axis.text.y = element_text(size = 10, vjust = 0.5, hjust=0.5))
    ggsave(paste0("../../plots/ccc/integrated/gojo_source-", val, "_targets-malignant_group.pdf"))
}


for (val in malignant_groups)
{
    print(val)
    gojo_ss2_liana %>%
    dplyr::filter(source ==val) %>%
    dplyr::top_n(25, desc(aggregate_rank)) %>%
    liana_dotplot(source_groups = c(val),
    target_groups = c_types)+
    theme(strip.text = element_text(size = 10, colour = "gray6"),
    plot.title = element_text(size=4), 
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust=0.5), 
    axis.text.y = element_text(size = 10, vjust = 0.5, hjust=0.5))
    ggsave(paste0("../../plots/ccc/integrated/gojo_source-", val, "_targets-mye_oli_stro_lymp.pdf"))
}
