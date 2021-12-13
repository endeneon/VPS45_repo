# Siwei 22 Jun 2021
# Reanalyse the CROP-Seq NGN2-Glut data by using cell clusters
# use NGN2_neuron_agg/R/NGN2_source.RData

library(Seurat)
library(sctransform)
library(glmGamPoi)


### There might be conflicts between packages,
### if running SCTransform, only load the three libraries above
##############################################################
library(ggplot2)
library(patchwork)


library(data.table)
library(ggplot2)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)

library(future)

library(readr)

###
## !! save the results to a different file!
load("NGN2_neuron_agg/R/NGN2_source.RData")

###
# split gRNA rows from Seurat matrix for subsequent quantification
NGN2_neuron_gRNA_stripped <- 
  subset(x = NGN2_neuron_with_gRNA, 
         features = rownames(NGN2_neuron_with_gRNA)[1:(nrow(NGN2_neuron_with_gRNA) - 76)])
# set another Seurat matrix contains gRNA counts only
NGN2_neuron_gRNA_only <- 
  subset(x = NGN2_neuron_with_gRNA, 
         features = rownames(NGN2_neuron_with_gRNA)[(nrow(NGN2_neuron_with_gRNA) - 75):
                                                      nrow(NGN2_neuron_with_gRNA)])

# make a data frame with each cell assigned to one unique gRNA identity
cell_gRNA_identity <- 
  rownames(NGN2_neuron_gRNA_only)[colMaxs(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))]
cell_gRNA_identity <- 
  as.data.frame(rbind(cell_gRNA_identity), 
                stringsAsFactors = F)
colnames(cell_gRNA_identity) <- 
  colnames(NGN2_neuron_gRNA_only)
# check number of cells assigned to each gRNA
table(rownames(NGN2_neuron_gRNA_only)[colMaxs(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))])

# perform dimensionality reduction, use PCA first then t-SNE, umap
NGN2_neuron_gRNA_stripped <- 
  FindVariableFeatures(NGN2_neuron_gRNA_stripped, 
                       verbose = T)
NGN2_neuron_gRNA_stripped <- RunPCA(NGN2_neuron_gRNA_stripped, 
                                    # features = 3000,
                                    verbose = T)
# check the dimensionality of signal, use PC 1-20
# ElbowPlot(NGN2_neuron_gRNA_stripped)
# tsne
NGN2_neuron_gRNA_stripped <- RunTSNE(NGN2_neuron_gRNA_stripped, 
                                     dims = 1:20, 
                                     verbose = T)
NGN2_neuron_gRNA_stripped <- RunUMAP(NGN2_neuron_gRNA_stripped, 
                                     dims = 1:20,
                                     verbose = T)
# find neighbours and cluster
NGN2_neuron_gRNA_stripped <- FindNeighbors(NGN2_neuron_gRNA_stripped,
                                           dims = 1:20,
                                           do.plot = T, 
                                           k.param = 5,
                                           nn.method = "rann",
                                           reduction = "pca",
                                           verbose = T)
NGN2_neuron_gRNA_stripped <- FindClusters(NGN2_neuron_gRNA_stripped, 
                                          resolution = 0.2)

# check plots
DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        label = T,
        repel = T)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("MAP2"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SLC17A6",
                         "SLC17A7"),
            blend = T,
            ncol = 2, 
            cols = c("red", "blue"))
# DimPlot(NGN2_neuron_gRNA_stripped, reduction = "tsne")
# check gene expression
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SLC17A6", "SLC17A7"), 
            blend = T)

StackedVlnPlot(obj = NGN2_neuron_gRNA_stripped,
               features = c("SLC17A6", "SLC17A7", 
                            "GAD2", "DLG4", "GLS",
                            "TH", 
                            "VIM", "SOX2")) +
  coord_flip()


save.image(file = "NGN2_reanalyse_no_MAP2.RData")

## assign gRNA identity to each cell
NGN2_neuron_gRNA_stripped$gRNA.identity <- 
  rownames(NGN2_neuron_gRNA_only)[colMaxs(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))]

## plot VPS452/3 gRNA
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA <- "unassigned"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "VPS45-1-gene"] <- "VPS45_gRNA"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "VPS45-2-gene"] <- "VPS45_gRNA"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "VPS45-3-gene"] <- "VPS45_gRNA"

NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-CTRL-00018-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-CTRL-00022-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-1-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-2-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-3-gene"] <- "ctrl"

NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA <-
  as.factor(NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA)

# set active identity
Idents(object = NGN2_neuron_gRNA_stripped) <- "plot_VPS45_2_3_gRNA"

DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        group.by = "plot_VPS45_2_3_gRNA",
        # cols = c("lightgrey", "blue", "red", "red"),
        order = c("unassigned", "ctrl",
                  "VPS45_gRNA"),
        cells.highlight = list(WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("ctrl")),
                               WhichCells(NGN2_neuron_gRNA_stripped,
                                          idents = c("VPS45_gRNA"))),
        # shuffle = T,
        sizes.highlight = .5,
        cols.highlight = c("red", "blue"))


DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        group.by = "plot_VPS45_2_3_gRNA",
        # cols = c("lightgrey", "blue", "red", "red"),
        order = c("unassigned", "ctrl",
                  "VPS45_gRNA"),
        cells.highlight = list(WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("ctrl")),
                               WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("VPS45_2_gRNA")),
                               WhichCells(NGN2_neuron_gRNA_stripped,
                                          idents = c("VPS45_3_gRNA"))),
        # shuffle = T,
        sizes.highlight = .5,
        cols.highlight = c("red", "red", "blue"))

##### Fig S3
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("GAD2"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("DLG4"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("GLS"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("TH"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("VIM"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SOX2"),
            blend = F)

#### plot cells assigned to unique/non-unique gRNAs
# make a data frame with each cell assigned to one unique gRNA identity
gRNA_stacked_plot_nonunique <- 
  rowSums(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))

# Select cells that have one gRNA expression more than 3x of all others
# 6658 cells passed this filter
NGN2_neuron_seurat_split <- NGN2_neuron_with_gRNA
NGN2_neuron_seurat_split <- NGN2_neuron_seurat_split[, 
                                         (colMaxs(as.matrix(NGN2_neuron_seurat_split@assays$RNA@counts[
                                           ((NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1] - 75):
                                              NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1]), ]), value = T) >
                                            Matrix::colSums(NGN2_neuron_seurat_split@assays$RNA@counts[
                                              ((NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1] - 75):
                                                 NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1]), ]) * 3 / 4)]
ncol(NGN2_neuron_seurat_split) 
NGN2_neuron_split_gRNA_only <- 
  subset(x = NGN2_neuron_seurat_split, 
         features = rownames(NGN2_neuron_seurat_split)[(nrow(NGN2_neuron_seurat_split) - 75):
                                                         nrow(NGN2_neuron_seurat_split)])

gRNA_stacked_plot_unique <- 
  rowSums(as.matrix(NGN2_neuron_split_gRNA_only@assays$RNA@counts))
