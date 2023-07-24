# Siwei 20 Dec 2022
# remove sample H12_sg_2 for enhanced DE genes
# Revision data for VPS45 CSC on rs2027349 AA/GG and VPS45/lncRNA/C1orf54
# synergistic effects

# 1. Compare enrichment pattern of AA/GG and AA/3xKD by WebGestaltR()
# 2. Check for synergistic effects of 1x/3x KD

# init
library(readr)

library(Rfast)
library(factoextra)
library(MASS)

pacman::p_load(limma, edgeR, pheatmap, RColorBrewer,
               ggplot2, ggpubr, qvalue, plyr, wesanderson,
               GSEABase, grid, scales, WebGestaltR, stringr)
# source functions
source(file = "Synergy_functions_only.R")

## Run DEG to get DE genes for 3xKD vs AA and GG vs AA
## need to test which GG/AA pair fits better

gene_count_table <-
  read_delim("VPS45/VPS45_3xKD_ReadsPerGene_STAR_08Dec2022_stranded_correct.txt",
             "\t",
             escape_double = FALSE,
             trim_ws = TRUE)

## Use ENSG gene id
df_raw <- gene_count_table
df_raw$Geneid <-
  str_split(string = df_raw$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_raw <-
  df_raw[!duplicated(df_raw$Geneid), ]
rownames(df_raw) <- df_raw$Geneid
####

# gencode_v35_ENSG_Genename_final <-
#   read_delim("gencode.v35.ENSG.Genename.final.list",
#              delim = "\t", escape_double = FALSE,
#              col_names = FALSE, trim_ws = TRUE)
#
# gencode_v35_ENSG_Genename_final <-
#   gencode_v35_ENSG_Genename_final[, c(1, 5)]
# colnames(gencode_v35_ENSG_Genename_final) <-
#   c("Geneid", "Gene_symbol")
#
# gencode_v35_ENSG_Genename_final <-
#   gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
# gencode_v35_ENSG_Genename_final <-
#   gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]
#
# df_raw <-
#   merge(x = gene_count_table,
#         y = gencode_v35_ENSG_Genename_final,
#         by = "Geneid")

# df_raw <-
#   df_raw[!duplicated(df_raw$Gene_symbol), ]
# rownames(df_raw) <- df_raw$Gene_symbol

df_4_DGE <- df_raw
df_4_DGE$Geneid <- NULL
# df_4_DGE$Gene_symbol <- NULL

rownames(df_4_DGE) <- df_raw$Geneid

df_4_DGE$`A11-EGFP-1_` <- NULL
df_4_DGE$`A11-EGFP-2_` <- NULL

df_4_DGE <-
  df_4_DGE[, c(4, 5, 6,
               7, 8, 9,
               1, 2, 3)]
df_4_DGE <-
  df_4_DGE[, -5]

rownames(df_4_DGE) <- df_raw$Geneid

# df_4_DGE$H12_EGF_1_ <- NULL # H12_EGF_1 showed as an outlier on MDS

## make metadata df
df_metadata <-
  data.frame(samples = colnames(df_4_DGE))
df_metadata$genotype <-
  factor(x = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
               "AA_3xKD", "AA_3xKD",
               "GG_EGFP", "GG_EGFP", "GG_EGFP"),
         levels = c("AA_EGFP", "AA_3xKD", "GG_EGFP"))
df_metadata$cell_line <-
  factor(x = c(rep_len("H12",
                       length.out = 5),
               rep_len("B11",
                       length.out = 3)),
         levels = c("H12", "B11"))

# df_metadata$genotype <-
#   factor(x = c("AA_EGFP", "AA_EGFP",
#                "AA_3xKD", "AA_3xKD", "AA_3xKD",
#                "GG_EGFP", "GG_EGFP", "GG_EGFP"),
#          levels = c("AA_EGFP", "AA_3xKD", "GG_EGFP"))

### make DGEList
df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = df_metadata$samples,
          group = df_metadata$genotype,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE),
                stringsAsFactors = F)
# hist(rowmeans(as.matrix(cpm_gene_count)),
#      breaks = 10000)


cpm_cutoff <- 1 # 17431 genes left

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
df_DGE <-
  df_DGE[(rowSums(cpm_gene_count[, df_metadata$genotype %in%
                                   "GG_EGFP"] >= cpm_cutoff) >= 2) |
           (rowSums(cpm_gene_count[, df_metadata$genotype %in%
                                     "AA_EGFP"] >= cpm_cutoff) >= 2) |
           (rowSums(cpm_gene_count[, df_metadata$genotype %in%
                                     "AA_3xKD"] >= cpm_cutoff) >= 2)
         , ]

retained_gene_list <- rownames(df_DGE)

df_DGE <- calcNormFactors(df_DGE)
plotMDS(df_DGE)

cpm_gene_count_filtered <-
  cpm_gene_count[rownames(cpm_gene_count) %in%
                   rownames(df_DGE)
                 , ]

df_design_matrix <-
  model.matrix(~ factor(genotype),
               data = df_metadata)

df_DGE <-
  estimateDisp(df_DGE,
               design = df_design_matrix)

results_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix)

results_QLM_3xKD_vs_AA <-
  glmQLFTest(results_QLM,
             coef = 2)
results_QLM_3xKD_vs_AA <-
  results_QLM_3xKD_vs_AA$table
results_QLM_3xKD_vs_AA$FDR <-
  p.adjust(results_QLM_3xKD_vs_AA$PValue,
           method = "fdr")

results_QLM_GG_vs_AA <-
  glmQLFTest(results_QLM,
             coef = 3)
results_QLM_GG_vs_AA <-
  results_QLM_GG_vs_AA$table
results_QLM_GG_vs_AA$FDR <-
  p.adjust(results_QLM_GG_vs_AA$PValue,
           method = "fdr")

save(list = c("results_QLM_3xKD_vs_AA",
              "results_QLM_GG_vs_AA"),
     file = "df_results_GG_3xKD_AA_no_sg2.RData")

