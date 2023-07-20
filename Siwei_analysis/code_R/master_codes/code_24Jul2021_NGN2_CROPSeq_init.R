# Siwei 21 Jul 2021
# Use SCTransform (note, install directly and compile locally by
# install.package() will work)

# Siwei 25 Oct 2019
# initialise dataset by Seurat
# analyse NGN2-neuron only (line 08, 09, and 11)
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

# set threads and parallelization
# this one will cause Cstack error if used with SCTransform() or NormalizeData()
plan("multicore", workers = 4)
plan("sequential")
plan()
options(expressions = 20000)

# plan("sequential")

# load data
NGN2_neuron_seurat <- Read10X(data.dir = "NGN2_neuron_agg/outs/filtered_feature_bc_matrix/")
NGN2_neuron_seurat <- CreateSeuratObject(counts = NGN2_neuron_seurat, 
                                      project = "NGN2_neuron_all")

nrow(NGN2_neuron_seurat) # 29922 genes
ncol(NGN2_neuron_seurat) # 30884 cells
rownames(NGN2_neuron_seurat)[(nrow(NGN2_neuron_seurat) - 75):nrow(NGN2_neuron_seurat)]
# 24404 cells with valid gRNAs
sum(Matrix::colSums(NGN2_neuron_seurat[(nrow(NGN2_neuron_seurat) - 75):nrow(NGN2_neuron_seurat), ]) > 0)

# data preprocessing
## remove cells without any gRNA reads
## 24404 cells with valid gRNAs
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


## Plot Fig. 2B from here, separate MAP2, SLC17A6, SLC17A7
# ## run SCTransform
# NGN2_neuron_with_gRNA <- SCTransform(NGN2_neuron_with_gRNA,
#                                      method = "glmGamPoi",
#                                      vars.to.regress = "percent.mt",
#                                      verbose = T)


# the SCTransform() function causes stack error, use the traditional approach
# including NormalizeData(), ScaleData(), and FindVariableFeatures
NGN2_neuron_with_gRNA <- NormalizeData(NGN2_neuron_with_gRNA, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = 10000,
                                       verbose = T)

# remove genes with all 0 count, 29922 -> 22351
NGN2_neuron_with_gRNA <- subset(x = NGN2_neuron_with_gRNA, 
                                features = rownames(NGN2_neuron_with_gRNA)[rowMaxs(as.matrix(NGN2_neuron_with_gRNA@assays$RNA@counts), value = T) > 0])

NGN2_neuron_with_gRNA <- ScaleData(NGN2_neuron_with_gRNA, 
                                   features = rownames(NGN2_neuron_with_gRNA),
                                   vars.to.regress = c("orig.ident",
                                                       "nCount_RNA", 
                                                       "percent.mt"), 
                                   verbose = T, 
                                   model.use = "negbinom")
# check gRNA efficiency
NGN2_neuron_with_gRNA <- FindVariableFeatures(NGN2_neuron_with_gRNA, verbose = T)
NGN2_neuron_with_gRNA <- RunPCA(NGN2_neuron_with_gRNA, verbose = T)
# ElbowPlot(NGN2_neuron_with_gRNA)
# NGN2_neuron_with_gRNA <- RunTSNE(NGN2_neuron_with_gRNA, 
#                                  tsne.method = "Rtsne")
# NGN2_neuron_with_gRNA <- RunUMAP(NGN2_neuron_with_gRNA, 
#                                  dims = 1:30,
#                                  verbose = T)



save.image(file = "NGN2_neuron_agg/R/NGN2_source.RData")
# find neighbours and clusters
## here use nn.method = "rann" or it may cause 
## "node stack overflow" error
NGN2_neuron_with_gRNA <- FindNeighbors(NGN2_neuron_with_gRNA,
                                       k.param = 5,
                                       nn.method = "rann",
                                       reduction = "pca",
                                       dims = 1:30)
## here FindCluster() requires single-thread
NGN2_neuron_with_gRNA <- FindClusters(NGN2_neuron_with_gRNA,
                                      resolution = 0.2, 
                                      verbose = T)
NGN2_neuron_with_gRNA <- RunTSNE(NGN2_neuron_with_gRNA, 
                                 tsne.method = "Rtsne")
NGN2_neuron_with_gRNA <- RunUMAP(NGN2_neuron_with_gRNA, 
                                 dims = 1:30,
                                 verbose = T)
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
