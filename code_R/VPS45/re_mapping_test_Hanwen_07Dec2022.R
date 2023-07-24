# Siwei 06 Dec 2022
# edgeR analysis for VPS45 3xKD RNASeq data

# init
library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)

library(ggplot2)
library(gplots)
library(RColorBrewer)


# load data
gene_count_table <- 
  read_delim("~/NVME/RNASeq_12Feb2021/STAR_output/Alena_GEO/output/ReadsPerGene_STAR.txt", 
             "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE)

df_raw <-
  gene_count_table

df_raw$Geneid <-
  str_split(df_raw$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]

df_raw <-
  df_raw[!duplicated(df_raw$Geneid), ]

df_4_DGE <- df_raw
df_4_DGE$Geneid
df_4_DGE$Geneid <- NULL
rownames(df_4_DGE) <- df_raw$Geneid

### make DGEList
df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = colnames(df_4_DGE),
          # group = df_metadata$genotype,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

DGE_cpm <-
  as.data.frame(cpm(df_4_DGE))

DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000285184", ] # lncRNA
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000136631", ] # VPS45
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000118292", ] # C1orf54

DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000284202", ] 
