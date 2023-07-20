# Siwei 16 Dec 2022
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
rownames(df_4_DGE) <- df_raw$Geneid

# df_4_DGE$H12_EGF_1_ <- NULL # H12_EGF_1 showed as an outlier on MDS

## make metadata df
df_metadata <-
  data.frame(samples = colnames(df_4_DGE))
df_metadata$genotype <-
  factor(x = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
               "AA_3xKD", "AA_3xKD", "AA_3xKD",
               "GG_EGFP", "GG_EGFP", "GG_EGFP"),
         levels = c("AA_EGFP", "AA_3xKD", "GG_EGFP"))
df_metadata$cell_line <-
  factor(x = c(rep_len("H12",
                       length.out = 6),
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
     file = "df_results_GG_3xKD_AA.RData")


### import logFC tables and merge
load("df_results_GG_3xKD_AA.RData")

singleKD_C1orf54_vs_AA <-
  read_delim("/nvmefs/VPS45/R_VPS45/Rapid_neuron_gene_exp_table_C1orf54KD_vs_EGFP_31Mar2021.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
singleKD_C1orf54_vs_AA <-
  singleKD_C1orf54_vs_AA[, c("Geneid", "logFC")]
singleKD_lncRNA_vs_AA <-
  read_delim("/nvmefs/VPS45/R_VPS45/Rapid_neuron_gene_exp_table_lncRNAKD_vs_EGFP_31Mar2021.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
singleKD_lncRNA_vs_AA <-
  singleKD_lncRNA_vs_AA[, c("Geneid", "logFC")]
singleKD_VPS45_vs_AA <-
  read_delim("/nvmefs/VPS45/R_VPS45/Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_31Mar2021.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
singleKD_VPS45_vs_AA <-
  singleKD_VPS45_vs_AA[, c("Geneid", "logFC")]

COS_RNAseq_SCZ_BD_GVEX <-
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/results/tables/RNAseq_SCZ_BD_GVEX.csv")
colnames(COS_RNAseq_SCZ_BD_GVEX)[1] <- "Geneid"
COS_RNAseq_SCZ <-
  COS_RNAseq_SCZ_BD_GVEX[, c(1, 2)]

CMC_RNAseq_SCZ_BD_CMC <-
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/results/tables/RNAseq_SCZ_BD_CMC.csv")
colnames(CMC_RNAseq_SCZ_BD_CMC)[1] <- "Geneid"
CMC_RNASeq_SCZ <-
  CMC_RNAseq_SCZ_BD_CMC[, c(1, 2)]
CMC_RNASeq_BIP <-
  CMC_RNAseq_SCZ_BD_CMC[, c(1, 10)]

RNAseq_ASD_4region_sumstats <-
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/results/tables/RNAseq_ASD_4region_sumstats.csv")
colnames(RNAseq_ASD_4region_sumstats)[1] <- "Geneid"
RNAseq_ASD_4region_sumstats <-
  RNAseq_ASD_4region_sumstats[, c(1, 2)]

Microarray_MDD_metaanalysis_092017 <-
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/results/tables/Microarray_MDD_metaanalysis_092017.csv")
colnames(Microarray_MDD_metaanalysis_092017)[1] <- "Geneid"
Microarray_MDD_metaanalysis_092017 <-
  Microarray_MDD_metaanalysis_092017[, c(1, 2)]

Microarray_AAD_metaanalysis_092017 <-
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/results/tables/Microarray_AAD_metaanalysis_092017.csv")
colnames(Microarray_AAD_metaanalysis_092017)[1] <- "Geneid"
Microarray_AAD_metaanalysis_092017 <-
  Microarray_AAD_metaanalysis_092017[, c(1, 2)]


# save.image(file = "all_logFC_data_sets.RData")
load("all_logFC_data_sets.RData")
load("df_results_GG_3xKD_AA.RData")

# merge data table
results_QLM_3xKD_vs_AA$Geneid <- rownames(results_QLM_3xKD_vs_AA)
results_QLM_3xKD_vs_AA <-
  results_QLM_3xKD_vs_AA[, c("Geneid", "logFC")]
results_QLM_GG_vs_AA$Geneid <- rownames(results_QLM_GG_vs_AA)
results_QLM_GG_vs_AA <-
  results_QLM_GG_vs_AA[, c("Geneid", "logFC")]

##
# build master data table
# remember to assign colnames at each step to differentiate logFC origin

df_master_logFC <-
  merge(x = results_QLM_3xKD_vs_AA,
        y = results_QLM_GG_vs_AA,
        by = "Geneid")
colnames(df_master_logFC) <-
  c("Geneid", "KDx3_vs_AA", "GG_vs_AA")
df_master_logFC <-
  merge(x = df_master_logFC,
        y = COS_RNAseq_SCZ,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "COS_RNAseq_SCZ"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = CMC_RNASeq_SCZ,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "CMC_RNASeq_SCZ"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = CMC_RNASeq_BIP,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "CMC_RNASeq_BIP"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = RNAseq_ASD_4region_sumstats,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "RNASeq_ASD"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = Microarray_MDD_metaanalysis_092017,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "Microarray_MDD"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = Microarray_AAD_metaanalysis_092017,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "Microarray_AAD"


df_master_logFC <-
  merge(x = df_master_logFC,
        y = singleKD_C1orf54_vs_AA,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "singleKD_C1orf54_vs_AA"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = singleKD_lncRNA_vs_AA,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "singleKD_lncRNA_vs_AA"

df_master_logFC <-
  merge(x = df_master_logFC,
        y = singleKD_VPS45_vs_AA,
        by = "Geneid")
colnames(df_master_logFC)[ncol(df_master_logFC)] <-
  "singleKD_VPS45_vs_AA"


save.image(file = "df_master_logFC.RData")

### split df_master_logFC into disease and calculated
df_master_logFC <-
  df_master_logFC[!duplicated(df_master_logFC$Geneid), ]

df_diseases <-
  df_master_logFC[, c(4, 5, 6, 7, 8, 9)]
rownames(df_diseases) <-
  df_master_logFC$Geneid
df_KD_results <-
  df_master_logFC[, c(10, 11, 12, 2, 3)]
rownames(df_KD_results) <-
  df_master_logFC$Geneid

########
## COS SCZ
cor.test(x = df_master_logFC$singleKD_C1orf54_vs_AA,
         y = df_master_logFC$COS_RNAseq_SCZ,
         method = "pearson",
         alternative = "t")

cor.test(x = df_master_logFC$KDx3_vs_AA,
         y = df_master_logFC$GG_vs_AA,
         method = "pearson",
         alternative = "t")

cor_return <-
  cor.test(x = df_master_logFC$singleKD_C1orf54_vs_AA,
           y = df_master_logFC$COS_RNAseq_SCZ,
           method = "pearson",
           alternative = "t")

calc_corr <-
  function(x, y) {
    pearson_cor_results <-
      cor.test(x = x,
               y = y,
               method = "pearson",
               alternative = "t")
    return_val <-
      vector(mode = "numeric",
             length = 2L)
    return_val[1] <- pearson_cor_results$estimate
    return_val[2] <- pearson_cor_results$p.value
    return(return_val)
  }

df_to_plot <-
  data.frame(Samples = c(1:ncol(df_KD_results)),
             Correlation = c(1:ncol(df_KD_results)),
             PVal = c(1:ncol(df_KD_results)))

for (i in 1:ncol(df_KD_results)) {

  corr_return <-
    cor.test(x = df_diseases[[1]],
             y = df_KD_results[[i]],
             method = "spearman",
             alternative = "t")
  # print(corr_return)
  df_to_plot$Samples[i] <- colnames(df_KD_results)[i]
  df_to_plot$Correlation[i] <- corr_return$estimate
  df_to_plot$PVal[i] <- corr_return$p.value
}
