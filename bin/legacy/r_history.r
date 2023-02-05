


if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
install.packages("remotes")
remotes::install_github('saezlab/liana')
library(tidyverse)
library(magrittr)
library(liana)
install.packages("tidyverse")
library(tidyverse)
library(magrittr)
library(liana)
show_resources()
show_methods()
setwd("Google Drive/Projects/saezlab/snRNASeq_Brain_Tumor_Data_Analysis/data/aniello_processed_objects/")
testdata <- readRDS(file ="human_Seurat_final_annotation.RData")
testdata <- load("human_Seurat_final_annotation.RData")
testdata %>% dplyr::glimpse()
testdata <- sce.m
testdata %>% dplyr::glimpse()
liana_test <- liana_wrap(testdata)
liana_test <- liana_wrap(testdata,idents_col = "final_annotation")
liana_test <- liana_test %>%
liana_aggregate()
dplyr::glimpse(liana_test)
liana_test %>%
dplyr::filter(source =="B") %>%
dplyr::top_n(25, desc(aggregate_rank))
sce.m$final_annotation
sce.m$final_annotation
liana_test$source
liana_test %>%
dplyr::filter(source =="B") %>%
dplyr::top_n(25, desc(aggregate_rank)) %>%
liana_dotplot(source_groups = c("malignant_B062_007+pdx01-x_human"),
target_groups = c("stromal", "myeloid T", "oligodendrocytes"))
liana_test %>%
dplyr::filter(source =="malignant_B062_007+pdx01-x_human") %>%
dplyr::top_n(25, desc(aggregate_rank)) %>%
liana_dotplot(source_groups = c("malignant_B062_007+pdx01-x_human"),
target_groups = c("stromal", "myeloid T", "oligodendrocytes"))
liana_test %>%
dplyr::filter(source =="malignant_B062_007+tumor01") %>%
dplyr::top_n(25, desc(aggregate_rank)) %>%
liana_dotplot(source_groups = c("malignant_B062_007+tumor01"),
target_groups = c("stromal", "myeloid T", "oligodendrocytes"))
liana_test %>%
dplyr::filter(source =="malignant_B062_007+tumor01") %>%
dplyr::top_n(25, desc(aggregate_rank)) %>%
liana_dotplot(source_groups = c("malignant_B062_007+tumor01"),
target_groups = c("stromal", "myeloid Q", "oligodendrocytes"))
liana_test
liana_test
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
history()
savehistory("~/Desktop/Untitled.Rhistory")
