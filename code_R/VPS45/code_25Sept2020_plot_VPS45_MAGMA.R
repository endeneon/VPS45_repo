# plot VPS45 MAGMA SNP enrichment
# Window @ [20k, 5k]
# Siwei 25 Sept 2020

# init
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(readxl)
library(colorRamps)

VPS45_GG_baseline_enrichment <- read_excel("/nvmefs/VPS45/MAGMA_20k_5k/table/VPS45_GG_baseline_enrichment.xlsx")

VPS45_GG_heatmat_matrix <- matrix(VPS45_GG_baseline_enrichment$`-logP`, 
                                  nrow = 3)
colnames(VPS45_GG_heatmat_matrix) <- unique(VPS45_GG_baseline_enrichment$Disease)
rownames(VPS45_GG_heatmat_matrix) <- VPS45_GG_baseline_enrichment$VARIABLE[1:3]

colnames(VPS45_GG_heatmap_matrix_v2)[5] <- "Crohn's Disease"

heatmap.2(as.matrix(VPS45_GG_heatmat_matrix),
          Rowv = F, Colv = F, 
          dendrogram = "none",
          scale = "none",
          col = matlab.like2(200),
          trace = "none",
          density.info = "none",
          margins = c(10,10),
          key.title = NA,
          key.xlab = "-log10P value",
          keysize = 2,
          srtCol = 45,
          cexRow = 1,
          cexCol = 1)


### v2

rownames(VPS45_GG_heatmap_matrix_v2) <- 
  rownames(VPS45_GG_heatmat_matrix)

heatmap.2(as.matrix(VPS45_GG_heatmap_matrix_v2),
          Rowv = F, Colv = F, 
          dendrogram = "none",
          scale = "none",
          col = matlab.like2(50),
          trace = "none",
          density.info = "none",
          margins = c(6,8),
          key.title = NA,
          key.xlab = "-log10P value",
          keysize = 2.3,
          srtCol = 45,
          cexRow = 1,
          cexCol = 1)
