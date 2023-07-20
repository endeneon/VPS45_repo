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


# library(data.table)
library(ggplot2)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)


library(readr)


# load("/nvmefs/VPS45/R_VPS45/NGN2_neuron_agg/R/results_per_gRNA.RData")


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
DimPlot(NGN2_neuron_MAP2_endogenous, 
        reduction = "umap",
        label = T, 
        repel = T)
