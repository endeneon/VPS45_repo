# Siwei 21 Jul 2021
# Seemingly SCTransform still does not work 
# for a new unknown reason
# start from a transformed RData file

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)


### There might be conflicts between packages,
### if running SCTransform, only load the three libraries above
##############################################################


library(data.table)
library(ggplot2)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)

library(future)

library(readr)

load("NGN2_neuron_agg/R/NGN2_source.RData",
     verbose = T)
rm(list = ls(pattern = "^edgeR_MAP2"))
save.image(file = "NGN2_seurat_start.RData")

###
