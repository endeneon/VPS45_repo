# Siwei 01 Dec 2020
# Use 09+11 samples only
# count matrix regenerated using CellRanger 5.0
# extract cluster-specific cells for analysis
# do not use SCTransform()

# Use Nov30.RData
# make DeSeq2-based analysis

# init
library(Seurat)
library(data.table)
# library(ggplot2)
# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)
library(future)
library(edgeR)
library(DESeq2)

# set multithreads
library("BiocParallel")
register(MulticoreParam(12))

# load data
load("~/NVME/CROPSeq_neuron_new_analysis_23Nov2020/R_0911_new_analysis/Nov30.RData")
# remove downstream analysis outputs
rm(list = c("Glut_neuron", "gRNA_list", "i",
            "DGE_raw", "DGE_design", "DGE_glmFit", "DGE_glmFTest"))
setwd(".")

## analyse by gRNA
# assign gRNA identity by gene to each cell as meta.data
# merge neg-CTRL and neg-EGFP as one "neg_all_ctrl"
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- NULL
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- 
  NGN2_neuron_gRNA_stripped@meta.data$gRNA.indiv
NGN2_neuron_gRNA_stripped$gRNA.bygRNA[str_detect(string = NGN2_neuron_gRNA_stripped$gRNA.bygRNA, 
                                                 pattern = "neg-")] <- "neg_all_ctrl"

### Glut group: # 0, 1, 2, 3, 7, 8, 9, 10, 13 (SLC17A6 + SLC17A7)
### vst transformed data are meant to be used for clustering only, 
### not for DE analysis, either for DESeq or EdgeR
Glut_neuron <- 
  NGN2_neuron_gRNA_stripped[, NGN2_neuron_gRNA_stripped$seurat_clusters %in% 
                              c(0, 1, 2, 3, 7, 8, 9, 10, 13)]

## construct DeSeq2 object
count_DeSeq2 <- as.matrix(Glut_neuron@assays$RNA@counts)
colData_DeSeq2 <- data.frame(gRNA = Glut_neuron$gRNA.bygRNA)
rownames(colData_DeSeq2) <- colnames(count_DeSeq2)

dds_DeSeq2 <- DESeqDataSetFromMatrix(countData = count_DeSeq2,
                                     colData = colData_DeSeq2,
                                     design = ~ gRNA)
# dds_DeSeq2 <- estimateSizeFactors(dds_DeSeq2,
#                                   type = "iterate")

## remove all unnecessary variables and run garbage cleaning
rm(list = c("NGN2_neuron_gRNA_stripped", "NGN2_neuron_seurat", 
            "NGN2_neuron_gRNA_only", "NGN2_neuron_0911", "Glut_neuron",
            "count_DeSeq2"))
gc(verbose = T)

dds_DeSeq2 <- DESeq(dds_DeSeq2,
                    # fitType = "local",
                    sfType = "iterate",
                    parallel = T,
                    BPPARAM = MulticoreParam(workers = 48,
                                             progressbar = T))
res_output <- results(dds_DeSeq2,
                    contrast = c("gRNA", "VPS45-1-gene", "neg_all_ctrl"),
                    parallel = T)

save.image(file = "DeSeq2_output.RData")
