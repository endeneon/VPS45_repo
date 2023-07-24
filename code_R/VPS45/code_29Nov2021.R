# Siwei 29 Nov 2021
# plot corr plot using the two gene sets

# init
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(colorspace)

library(readxl)
library(readr)
library(stringr)
library(factoextra)

library(MASS)
library(mgcv)

####
df_gene_list_raw <-
  read_excel("Data_Fig6F.xlsx",
             sheet = 1)

df_gene_list_neuron_diff <- df_gene_list_raw[1:77, ]
df_gene_list_synaptic <- df_gene_list_raw[78:nrow(df_gene_list_raw), ]  


###

# load data

df_main <- read_excel("VPS45_AA_AG_GG_all_gene_table.xlsx", 
                      sheet = "VPS45_AA_AG_GG_all_gene_table_w")

# VPS45KD, scaling factor 1.329
df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,6)]
# add scaling factor
# df_to_merge$logFC <- df_to_merge$logFC + log2(1.329)
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_VPS45_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_VPS45_KD")
rm(df_to_merge)

# lncRNAKD, scaling factor 1.043
df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_lncRNAKD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,6)]
# df_to_merge$logFC <- df_to_merge$logFC + log2(1.043)
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_lncRNA_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_lncRNA_KD")
rm(df_to_merge)


# C1orf54KD, scaling factor 1.433
df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_C1orf54KD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,6)]
# df_to_merge$logFC <- df_to_merge$logFC + log2(1.433)
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_C1orf54_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_C1orf54_KD")
rm(df_to_merge)

df_main <- df_main[!duplicated(df_main$Geneid), ]


##### test simple additive model #####

# use GO term enriched list
raw_df <- df_main

# raw_df <- raw_df[raw_df$FDR < 0.05, ]
raw_df <- raw_df[raw_df$Gene_Symbol %in% df_gene_list_neuron_diff$Gene_Symbol, ]
raw_df <- raw_df[!duplicated(raw_df$Gene_Symbol), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

# use GO term enriched list
raw_df <- df_main

# raw_df <- raw_df[raw_df$FDR < 0.05, ]
raw_df <- raw_df[raw_df$Gene_Symbol %in% df_gene_list_synaptic$Gene_Symbol, ]
raw_df <- raw_df[!duplicated(raw_df$Gene_Symbol), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

# use GO term enriched list
raw_df <- df_main

# raw_df <- raw_df[raw_df$FDR < 0.05, ]
raw_df <- raw_df[raw_df$Gene_Symbol %in% df_gene_list_raw$Gene_Symbol, ]
raw_df <- raw_df[!duplicated(raw_df$Gene_Symbol), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

