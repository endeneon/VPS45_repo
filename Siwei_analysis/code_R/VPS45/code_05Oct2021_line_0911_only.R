# Siwei 05 Oct 2021
# !! only analyse line 09+11 !!!


# Siwei 21 Jul 2021
# Use SCTransform (note, install directly and compile locally by
# install.package() will work

# Siwei 25 Oct 2019
# initialise dataset by Seurat
# analyse NGN2-neuron only (09, and 11)
# perform further analysis by Seurat
# Since the original dataset were made using a filtered set of bc_matrix
# try not to impose further filters such as MT, nFeatures, etc.
# find the SLC17A6, SLC17A7, and MAP2 positive cells and 
# use those cells only for further analysis

# init
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(future)
# set threads and parallelization
# this one will cause Cstack error if used with SCTransform() or NormalizeData()
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 500000)

### There might be conflicts between packages,
### if running SCTransform, only load the three libraries above
### SCTransform is not compatible with Rfast
##############################################################
library(ggplot2)
library(patchwork)

library(RColorBrewer)
library(colorRamps)
library(scales)

# library(data.table)
# library(ggplot2)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)


library(readr)



# plan("sequential")

# load data
NGN2_neuron_seurat <- 
  Read10X(data.dir = "../../CROPSeq_neuron_new_analysis_23Nov2020/neuron_0911_new_23Nov2020/outs/filtered_feature_bc_matrix")
NGN2_neuron_seurat <- CreateSeuratObject(counts = NGN2_neuron_seurat, 
                                         project = "NGN2_neuron_all")

nrow(NGN2_neuron_seurat) # 29922 genes
ncol(NGN2_neuron_seurat) # 13270 cells
rownames(NGN2_neuron_seurat)[(nrow(NGN2_neuron_seurat) - 75):nrow(NGN2_neuron_seurat)]
# 10247 cells with valid gRNAs
sum(Matrix::colSums(NGN2_neuron_seurat[(nrow(NGN2_neuron_seurat) - 75):nrow(NGN2_neuron_seurat), ]) > 0)

# data preprocessing
## remove cells without any gRNA reads
## 10247 cells with valid gRNAs
NGN2_neuron_with_gRNA <- 
  NGN2_neuron_seurat[, 
                     Matrix::colSums(NGN2_neuron_seurat[(nrow(NGN2_neuron_seurat) - 75)
                                                        : nrow(NGN2_neuron_seurat), ]) > 0]
# Select cells that have one gRNA expression more than 3x of all others
# 8248 cells passed this filter
NGN2_neuron_with_gRNA <- 
  NGN2_neuron_with_gRNA[, 
                        (colMaxs(as.matrix(NGN2_neuron_with_gRNA@assays$RNA@counts[
                          ((NGN2_neuron_with_gRNA@assays$RNA@counts@Dim[1] - 75):
                             NGN2_neuron_with_gRNA@assays$RNA@counts@Dim[1]), ]), value = T) >
                           Matrix::colSums(NGN2_neuron_with_gRNA@assays$RNA@counts[
                             ((NGN2_neuron_with_gRNA@assays$RNA@counts@Dim[1] - 75):
                                NGN2_neuron_with_gRNA@assays$RNA@counts@Dim[1]), ]) * 3 / 4)]

# set additional properties
# cut out the last digit of each column name for cell identity (batch)
NGN2_neuron_with_gRNA@meta.data$orig.ident <- 
  as.factor(str_sub(colnames(NGN2_neuron_with_gRNA), -1))
# add a meta.data slot for MT gene percentage
NGN2_neuron_with_gRNA <- 
  PercentageFeatureSet(NGN2_neuron_with_gRNA, 
                       pattern = "^MT-", 
                       col.name = "percent.mt")
# detach "Rfast" before running SCTransform
detach("package:Rfast", unload = TRUE)

save.image(file = "NGN2_0911.RData")
load("/nvmefs/VPS45/R_VPS45/NGN2_0911.RData")
## Plot Fig. 2B from here, separate MAP2, SLC17A6, SLC17A7
## run SCTransform
NGN2_neuron_with_gRNA <- SCTransform(NGN2_neuron_with_gRNA,
                                     method = "glmGamPoi",
                                     vars.to.regress = "percent.mt",
                                     verbose = T)
# reattach "RF

# the SCTransform() function causes stack error, use the traditional approach
# including NormalizeData(), ScaleData(), and FindVariableFeatures
# NGN2_neuron_with_gRNA <- NormalizeData(NGN2_neuron_with_gRNA, 
#                                        normalization.method = "LogNormalize", 
#                                        scale.factor = 10000,
#                                        verbose = T)

# remove genes with all 0 count, 29922 -> 22351
NGN2_neuron_with_gRNA <- subset(x = NGN2_neuron_with_gRNA, 
                                features = rownames(NGN2_neuron_with_gRNA)[rowMaxs(as.matrix(NGN2_neuron_with_gRNA@assays$RNA@counts), value = T) > 0])

# NGN2_neuron_with_gRNA <- ScaleData(NGN2_neuron_with_gRNA, 
#                                    features = rownames(NGN2_neuron_with_gRNA),
#                                    vars.to.regress = c("orig.ident",
#                                                        "nCount_RNA", 
#                                                        "percent.mt"), 
#                                    verbose = T, 
#                                    model.use = "negbinom")
# check gRNA efficiency
NGN2_neuron_with_gRNA <- FindVariableFeatures(NGN2_neuron_with_gRNA, verbose = T)
NGN2_neuron_with_gRNA <- RunPCA(NGN2_neuron_with_gRNA, verbose = T)
# ElbowPlot(NGN2_neuron_with_gRNA)
# NGN2_neuron_with_gRNA <- RunTSNE(NGN2_neuron_with_gRNA, 
#                                  tsne.method = "Rtsne")
# NGN2_neuron_with_gRNA <- RunUMAP(NGN2_neuron_with_gRNA, 
#                                  dims = 1:30,
#                                  verbose = T)




# find neighbours and clusters
## here use nn.method = "rann" or it may cause 
## "node stack overflow" error
NGN2_neuron_with_gRNA <- RunTSNE(NGN2_neuron_with_gRNA, 
                                 tsne.method = "Rtsne")
NGN2_neuron_with_gRNA <- RunUMAP(NGN2_neuron_with_gRNA, 
                                 dims = 1:30,
                                 verbose = T)
NGN2_neuron_with_gRNA <- FindNeighbors(NGN2_neuron_with_gRNA,
                                       k.param = 5,
                                       nn.method = "rann",
                                       reduction = "pca",
                                       dims = 1:30)
## here FindCluster() requires single-thread
NGN2_neuron_with_gRNA <- FindClusters(NGN2_neuron_with_gRNA,
                                      resolution = 0.2, 
                                      verbose = T)


save.image(file = "NGN2_0911.RData")


DimPlot(NGN2_neuron_with_gRNA, 
        reduction = "umap",
        label = T, 
        repel = T)

FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("MAP2", 
                         "SLC17A6",
                         "SLC17A7"),
            blend = F)
FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("MAP2"),
            blend = F)
FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("SLC17A6"),
            blend = F)
FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("SLC17A7"),
            blend = F)
FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("SLC17A6",
                         "SLC17A7"),
            blend = T,
            ncol = 2, 
            cols = c("magenta", "blue"))

### make stacked violin plot
# recapitulate the scanpy stacked violin plot function in R
# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
### ! NOTE ! ggplot2 changed behaviour
### need to change           
# plot.title = element_blank(),
# axis.title.x = element_blank(),
### to avoid excessive white space
### main functions #####
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj,
                           feature,
                           pt.size = 0,
                           # plot.margin = unit(c(0,0,0,0), "pt"),
                           plot.margin = unit(c(-3, 0, -3, 0), "pt"),
                           ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(0.8), hjust = rel(0.1),
                                      vjust = rel(0.5), angle = 0),
          plot.margin = plot.margin,
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = rel(1)))
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot <- function(obj, features,
                           pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x +
                             scale_y_continuous(breaks = c(y)) +
                             expand_limits(y = y))# +
  # scale_x_discrete(limits = c("1", "3", "8", "9", "12",
  #                             "0", "4", "6", "10", "11",
  #                             "2", "5", "13",
  #                             "7")))
  
  p <- patchwork::wrap_plots(plotlist = plot_list, 
                             ncol = 1,
                             heights = 20,
                             widths = 20)
  return(p)
}

StackedVlnPlot(obj = NGN2_neuron_with_gRNA,
               features = c("SLC17A6", "SLC17A7", 
                            "GAD2", "DLG4", "GLS",
                            "TH", 
                            "VIM", "SOX2")) +
  coord_flip()



# select cells that are either SLC17A6 or SLC17A7-positive
# 3882 cells passed this filter
NGN2_neuron_with_gRNA <- NGN2_neuron_with_gRNA[, ((NGN2_neuron_with_gRNA@assays$RNA@counts[rownames(NGN2_neuron_with_gRNA@assays$RNA@counts) %in% "SLC17A6", ] > 0) |
                                                    (NGN2_neuron_with_gRNA@assays$RNA@counts[rownames(NGN2_neuron_with_gRNA@assays$RNA@counts) %in% "SLC17A7", ] > 0))
]


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


Idents(object = NGN2_neuron_gRNA_stripped) <- "orig.ident"
VlnPlot(NGN2_neuron_gRNA_stripped,
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        ncol = 3,
        pt.size = 0)

### expressions
# str_remove_all(Glut_subset$gRNA.identity,
#                pattern = "-[[:digit:]]-gene")

unique(Glut_subset$gRNA.bygene)
unique(Glut_subset$gRNA.indiv)

# BAG5
df_test <- 
  data.frame(category.name = Glut_subset$gRNA.bygene[Glut_subset$gRNA.bygene %in% 
                                                               c("neg_all_ctrl", "BAG5")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$gRNA.bygene %in% 
                                                     c("neg_all_ctrl", "BAG5")],
             stringsAsFactors = F)

df_test$read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5", 
                                 Glut_subset$gRNA.bygene %in% c("neg_all_ctrl", "BAG5")])

df_test$gRNA.name <- as.factor(df_test$gRNA.name)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("neg_all_ctrl",
                                           "BAG5"))


ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  xlab("") +
  ylab("BAG5 Normalized expression value") +
  scale_fill_manual(name = "", 
                    values = c("lightblue", "orchid")) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic() +
  theme(legend.position = "none")

# p-value = 0.2011
t.test(x = df_test$read.value[df_test$category.name %in% "neg_all_ctrl"], 
       y = df_test$read.value[df_test$category.name %in% "BAG5"], 
       var.equal = F,
       alternative = "t",
       paired = F)

# PBRM1
df_test <- 
  data.frame(category.name = Glut_subset$gRNA.bygene[Glut_subset$gRNA.bygene %in% 
                                                       c("neg_all_ctrl", "PBRM1")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$gRNA.bygene %in% 
                                                     c("neg_all_ctrl", "PBRM1")],
             stringsAsFactors = F)

df_test$read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "PBRM1", 
                                 Glut_subset$gRNA.bygene %in% c("neg_all_ctrl", "PBRM1")])

df_test$gRNA.name <- as.factor(df_test$gRNA.name)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("neg_all_ctrl",
                                           "PBRM1"))


ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  xlab("") +
  ylab("PBRM1 normalized expression value") +
  scale_fill_manual(name = "", 
                    values = c("lightblue", "orchid")) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic() +
  theme(legend.position = "none")

# p-value = 0.00148
t.test(x = df_test$read.value[df_test$category.name %in% "neg_all_ctrl"], 
       y = df_test$read.value[df_test$category.name %in% "PBRM1"], 
       var.equal = F,
       alternative = "t",
       paired = F)

# PBRM1 use PBRM1 gRNA #2
df_test <- 
  data.frame(category.name = Glut_subset$gRNA.bygRNA[Glut_subset$gRNA.bygRNA %in% 
                                                       c("neg_all_ctrl", "PBRM1-2-gene")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$gRNA.bygRNA %in% 
                                                     c("neg_all_ctrl", "PBRM1-2-gene")],
             stringsAsFactors = F)

df_test$read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "PBRM1", 
                                 Glut_subset$gRNA.bygRNA %in% c("neg_all_ctrl", "PBRM1-2-gene")])

df_test$gRNA.name <- as.factor(df_test$gRNA.name)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("neg_all_ctrl",
                                           "PBRM1-2-gene"))


ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  xlab("") +
  ylab("GNL3 normalized expression value") +
  scale_fill_manual(name = "", 
                    values = c("lightblue", "orchid")) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic() +
  theme(legend.position = "none")

# p-value = 0.01343
t.test(x = df_test$read.value[df_test$category.name %in% "neg_all_ctrl"], 
       y = df_test$read.value[df_test$category.name %in% "PBRM1-2-gene"], 
       var.equal = F,
       alternative = "t",
       paired = F)


# GNL3
df_test <- 
  data.frame(category.name = Glut_subset$gRNA.bygene[Glut_subset$gRNA.bygene %in% 
                                                       c("neg_all_ctrl", "PBRM1")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$gRNA.bygene %in% 
                                                     c("neg_all_ctrl", "PBRM1")],
             stringsAsFactors = F)

df_test$read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "GNL3", 
                                 Glut_subset$gRNA.bygene %in% c("neg_all_ctrl", "PBRM1")])

df_test$gRNA.name <- as.factor(df_test$gRNA.name)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("neg_all_ctrl",
                                           "PBRM1"))


ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  xlab("") +
  ylab("GNL3 normalized expression value") +
  scale_fill_manual(name = "", 
                    values = c("lightblue", "orchid")) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic() +
  theme(legend.position = "none")

# p-value = 0.00148
t.test(x = df_test$read.value[df_test$category.name %in% "neg_all_ctrl"], 
       y = df_test$read.value[df_test$category.name %in% "PBRM1"], 
       var.equal = F,
       alternative = "t",
       paired = F)

# GNL3 use PBRM1 gRNA #2
df_test <- 
  data.frame(category.name = Glut_subset$gRNA.bygRNA[Glut_subset$gRNA.bygRNA %in% 
                                                       c("neg_all_ctrl", "PBRM1-2-gene")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$gRNA.bygRNA %in% 
                                                     c("neg_all_ctrl", "PBRM1-2-gene")],
             stringsAsFactors = F)

df_test$read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "GNL3", 
                                 Glut_subset$gRNA.bygRNA %in% c("neg_all_ctrl", "PBRM1-2-gene")])

df_test$gRNA.name <- as.factor(df_test$gRNA.name)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("neg_all_ctrl",
                                           "PBRM1-2-gene"))


ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  xlab("") +
  ylab("GNL3 normalized expression value") +
  scale_fill_manual(name = "", 
                    values = c("lightblue", "orchid")) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic() +
  theme(legend.position = "none")

# p-value = 0.01504
t.test(x = df_test$read.value[df_test$category.name %in% "neg_all_ctrl"], 
       y = df_test$read.value[df_test$category.name %in% "PBRM1-2-gene"], 
       var.equal = F,
       alternative = "t",
       paired = F)


######

# make bar plot
load("/nvmefs/VPS45/R_VPS45/NGN2_0911.RData")
# set another Seurat matrix contains gRNA counts only
NGN2_neuron_nonunique <- 
  subset(x = NGN2_neuron_seurat, 
         features = rownames(NGN2_neuron_seurat)[(nrow(NGN2_neuron_seurat) - 75):
                                                   nrow(NGN2_neuron_seurat)])

NGN2_neuron_unique <- 
  subset(x = NGN2_neuron_with_gRNA, 
         features = rownames(NGN2_neuron_with_gRNA)[(nrow(NGN2_neuron_with_gRNA) - 75):
                                                      nrow(NGN2_neuron_with_gRNA)])


gRNA_stacked_plot <- 
  data.frame(gRNA = rownames(NGN2_neuron_nonunique@assays$RNA@counts),
             count_nonunique = rowSums(NGN2_neuron_nonunique@assays$RNA@counts),
             count_unique = rowSums(NGN2_neuron_unique@assays$RNA@counts),
             stringsAsFactors = T)

gRNA_stacked_plot$gRNA <-
  str_replace_all(gRNA_stacked_plot$gRNA,
                  pattern = "-gene",
                  replacement = "-gRNA")

ggplot(data = gRNA_stacked_plot, aes(x = gRNA)) +
  geom_bar(aes(y = gRNA_stacked_plot$count_nonunique), 
           stat = "identity", 
           colour = "#FFFFFF", 
           fill = "#888888") +
  geom_bar(aes(y = gRNA_stacked_plot$count_unique), 
           stat = "identity", 
           colour = "#FFFFFF", 
           fill = "#0000FF") +
  xlab("") +
  ylab("Cell count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 6, 
                                   colour = "black"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))
