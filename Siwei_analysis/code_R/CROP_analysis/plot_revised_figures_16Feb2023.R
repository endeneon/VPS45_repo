# make additional plots for CSC revision
# Siwei 16 Feb 2023


# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
plan("multisession", workers = 12)
# plan("sequential")
# plan()
options(expressions = 500000)
options(future.globals.maxSize = 2097152000)
### There might be conflicts between packages,
### if running SCTransform, only load the three libraries above
### SCTransform is not compatible with Rfast
##############################################################


# library(data.table)
library(ggplot2)
library(RColorBrewer)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)

library(edgeR)

library(readr)

# load data
load("Nov30.RData")
# remove downstream analysis outputs
rm(list = c("Glut_neuron", "gRNA_list", "i",
            "DGE_raw", "DGE_design", "DGE_glmFit", "DGE_glmFTest"))


# assign gRNA identity by gene to each cell as meta.data
# merge neg-CTRL and neg-EGFP as one "neg_all_ctrl"
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- NULL
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- 
  NGN2_neuron_gRNA_stripped@meta.data$gRNA.indiv
NGN2_neuron_gRNA_stripped$gRNA.bygRNA[str_detect(string = NGN2_neuron_gRNA_stripped$gRNA.bygRNA, 
                                                 pattern = "neg-")] <- "neg_all_ctrl"

DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        label = T, 
        repel = T)

CD_09_barcodes <- 
  read_csv("CD_09_barcodes.txt", 
           col_names = FALSE)
CD_09_barcodes <-
  CD_09_barcodes$X1

CD_11_barcodes <- 
  read_csv("CD_11_barcodes.txt", 
           col_names = FALSE)
CD_11_barcodes <-
  CD_11_barcodes$X1
length(CD_11_barcodes)

NGN2_neuron_gRNA_stripped@meta.data$line_ident <- ""
NGN2_neuron_gRNA_stripped@meta.data$line_ident[colnames(NGN2_neuron_gRNA_stripped) %in% CD_09_barcodes] <- "CD-09"
NGN2_neuron_gRNA_stripped@meta.data$line_ident[colnames(NGN2_neuron_gRNA_stripped) %in% CD_11_barcodes] <- "CD-11"

Idents(NGN2_neuron_gRNA_stripped) <- "line_ident"

DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        group.by = "line_ident",
        pt.size = 0.2,
        # cols = c("lightgrey", "blue", "red", "red"),
        # order = c("CD-09", "CD-11"),
        cells.highlight = list(WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("CD-09")),
                               WhichCells(NGN2_neuron_gRNA_stripped,
                                          idents = c("CD-11"))),
        # shuffle = T,
        sizes.highlight = .2,
        cols.highlight = c("blue", "red")) +
  scale_colour_manual(values = c("grey", "red", "blue"),
                      labels = c("Unidentified", "CD-11", "CD-09")) +
  ggtitle("Cells by line")

DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        group.by = "line_ident",
        pt.size = 0.2,
        # cols = c("lightgrey", "blue", "red", "red"),
        # order = c("CD-09", "CD-11"),
        cells.highlight = list(WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("CD-09"))),
        # shuffle = T,
        sizes.highlight = .2,
        cols.highlight = c("red", "blue")) +
  # scale_fill_discrete(labels = c("Unidentified", "CD-09", "CD-11")) +
  ggtitle("Grey=unidentified; Red = CD-09; Blue = CD-11")


### violin plots
# load data
load("Nov30.RData")
# remove previous downstream analysis outputs
rm(list = c("Glut_neuron", "gRNA_list", "i",
            "DGE_raw", "DGE_design", "DGE_glmFit", "DGE_glmFTest"))


# assign gRNA identity by gene to each cell as meta.data
# merge neg-CTRL and neg-EGFP as one "neg_all_ctrl"
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- NULL
NGN2_neuron_gRNA_stripped <- NGN2_neuron_0911

NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- 
  NGN2_neuron_gRNA_stripped@meta.data$gRNA.indiv

NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA[str_detect(string = NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA, 
                                                           pattern = "neg-")] <- "neg_all_ctrl"
CD_09_barcodes <- 
  read_csv("CD_09_barcodes.txt", 
           col_names = FALSE)
CD_09_barcodes <-
  CD_09_barcodes$X1

CD_11_barcodes <- 
  read_csv("CD_11_barcodes.txt", 
           col_names = FALSE)
CD_11_barcodes <-
  CD_11_barcodes$X1
length(CD_11_barcodes)

NGN2_neuron_gRNA_stripped@meta.data$line_ident <- ""
NGN2_neuron_gRNA_stripped@meta.data$line_ident[colnames(NGN2_neuron_gRNA_stripped) %in% CD_09_barcodes] <- "CD-09"
NGN2_neuron_gRNA_stripped@meta.data$line_ident[colnames(NGN2_neuron_gRNA_stripped) %in% CD_11_barcodes] <- "CD-11"

Idents(NGN2_neuron_gRNA_stripped) <- "line_ident"

### plotting
NGN2_neuron_with_gRNA <- SCTransform(NGN2_neuron_gRNA_stripped,
                                     method = "glmGamPoi",
                                     vars.to.regress = "percent.mt",
                                     variable.features.n = 10000, 
                                     verbose = T)
NGN2_neuron_with_gRNA <- RunTSNE(NGN2_neuron_with_gRNA, 
                                 tsne.method = "Rtsne", 
                                 seed.use = 42)
NGN2_neuron_with_gRNA <- RunUMAP(NGN2_neuron_with_gRNA, 
                                 reduction = "pca",
                                 dims = 1:30,
                                 verbose = T, 
                                 seed.use = 42)
NGN2_neuron_with_gRNA <- FindNeighbors(NGN2_neuron_with_gRNA,
                                       k.param = 5,
                                       nn.method = "rann",
                                       reduction = "pca",
                                       dims = 1:30)
## here FindCluster() requires single-thread
NGN2_neuron_with_gRNA <- FindClusters(NGN2_neuron_with_gRNA,
                                      resolution = 0.18, 
                                      graph.name = "RNA_snn",
                                      verbose = T)
DimPlot(NGN2_neuron_with_gRNA, 
        reduction = "umap",
        label = T, 
        repel = T)

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
