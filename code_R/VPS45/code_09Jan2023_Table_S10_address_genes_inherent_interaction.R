# Siwei 09 Jan 2023
# Address the reviewer's question:
# "I was wondering whether the conditioning (choice of genes put forward into 
# Table S10) could induce an interaction.The criterion for entry into Table S10 
# are genes found to significantly change expression of GG from CRISPR/Cas9-edited 
# NGN2-Glut (genotype AA vs. GG) which are derived from independent experiments to 
# the DEG derived from knockdown of VPS45, lncRNA and C1orf54 in A/A NGN2-Glut. I
# checked using simulation using the variance-covariance matrix estimated from
#Table S10, but maybe the authors can investigate this using realistic
# simulation as they data are new to me.

# init

library(readr)
library(readxl)

library(MASS)
library(mgcv)
library(lmtest)

library(ggplot2)
library(RColorBrewer)

# load data
load(file = "../R_synergy_analysis/VPS45_correlation_model_raw.RData")

# raw_df contains the 1267 most DE genes used in original Table S10
# df_main contains all genes that can be used

###

df_main <- read_excel("VPS45_AA_AG_GG_all_gene_table.xlsx", 
                      sheet = "VPS45_AA_AG_GG_all_gene_table_w")

# VPS45KD, scaling factor 1.329
df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,3,6)]
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
df_to_merge <- df_to_merge[, c(1,2,3,6)]
# df_to_merge$logFC <- df_to_merge$logFC + log2(1.043)
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_lncRNA_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_lncRNA_KD")
rm(df_to_merge)


# C1orf54KD, scaling factor 1.433
df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_C1orf54KD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,3,6)]
# df_to_merge$logFC <- df_to_merge$logFC + log2(1.433)
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_C1orf54_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_C1orf54_KD")
rm(df_to_merge)

df_main <- df_main[!duplicated(df_main$Geneid), ]
df_main_backup <- df_main

###

df_main <- 
  df_main[df_main$Geneid %in% df_test$Geneid, ]

df_test <-
  df_main[order(abs(df_main$logFC.groupGG),
                decreasing = T), ]

df_test <-
  df_test[1:1267, ]

# weighted_fit_test_df <-
#   lm(logFC.groupGG ~ 
#        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
#        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
#        (logFC_VPS45_KD * logFC_lncRNA_KD) +
#        (logFC_VPS45_KD * logFC_C1orf54_KD), 
#      data = df_test)

weighted_fit_test_df <-
  gam(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD), 
     data = df_test)

weighted_fit_test_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD), 
      data = raw_df)

weighted_fit_test_df <-
  lm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD), 
      data = raw_df)

coeftest(weighted_fit_test_df)

summary(weighted_fit_test_df)

### extract the non-interacting linear part
non_interacting_weighted_fit <-
  lm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
     data = df_main)

non_interacting_weighted_fit <-
  gam(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
     data = df_main)

coeftest(non_interacting_weighted_fit)
summary(non_interacting_weighted_fit)


#### make simple gam model use all genes
weighted_fit_test_df_single_factor <-
  gam(logFC.groupGG ~ 0 +
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD, 
      data = df_main_backup)
weighted_fit_test_df_single_factor <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD, 
      data = raw_df)

weighted_fit_test_df_single_factor <-
  lm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD), 
      data = raw_df)

weighted_fit_test_df_single_factor <-
  lm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD, 
     data = raw_df)
# coeftest(weighted_fit_test_df_single_factor)$Estimate
# weighted_fit_test_df_single_factor <- df_main_backup

y_prime <-
  raw_df

weighted_fit_test_df_single_factor <-
  as.data.frame(cbind(weighted_fit_test_df_single_factor$model,
                      weighted_fit_test_df_single_factor$residuals))
rownames(weighted_fit_test_df_single_factor) <- raw_df$Geneid
weighted_fit_test_df_single_factor <-
  weighted_fit_test_df_single_factor[rownames(weighted_fit_test_df_single_factor) %in%
                                       raw_df$Geneid, ]
weighted_fit_test_df_single_factor$y <-
  (weighted_fit_test_df_single_factor$coefficients[2] * raw_df$logFC_VPS45_KD) +
  (weighted_fit_test_df_single_factor$coefficients[4] * raw_df$logFC_C1orf54_KD) +
  (weighted_fit_test_df_single_factor$coefficients[3] * raw_df$logFC_lncRNA_KD)

# +
  weighted_fit_test_df_single_factor$`weighted_fit_test_df_single_factor$residuals`

weighted_fit_test_df_single_factor$y <-
  (weighted_fit_test_df_single_factor$logFC_VPS45_KD * raw_df$logFC_VPS45_KD) +
  (weighted_fit_test_df_single_factor$logFC_C1orf54_KD * raw_df$logFC_C1orf54_KD) +
  (weighted_fit_test_df_single_factor$logFC_lncRNA_KD * raw_df$logFC_lncRNA_KD) 

  
weighted_fit_test_real <-
  gam(y ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = weighted_fit_test_df_single_factor)

weighted_fit_test_real <-
  lm(y ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD), 
      data = weighted_fit_test_df_single_factor)
summary(weighted_fit_test_real)

weighted_fit_test_real <-
  lm(y ~ 
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD) -
       (logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD), 
     data = weighted_fit_test_df_single_factor)
summary(weighted_fit_test_real)



#### make simple gam model use all genes
weighted_fit_test_df_single_factor <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD, 
      data = df_main_backup)
# weighted_fit_test_df_single_factor <- df_main_backup

weighted_fit_test_df_single_factor <-
  as.data.frame(cbind(weighted_fit_test_df_single_factor$model,
                      weighted_fit_test_df_single_factor$residuals))
rownames(weighted_fit_test_df_single_factor) <- df_main_backup$Geneid
weighted_fit_test_df_single_factor <-
  weighted_fit_test_df_single_factor[rownames(weighted_fit_test_df_single_factor) %in%
                                       raw_df$Geneid, ]
weighted_fit_test_df_single_factor$y <-
  (weighted_fit_test_df_single_factor$logFC_VPS45_KD * df_main$logFC_VPS45_KD) +
  (weighted_fit_test_df_single_factor$logFC_C1orf54_KD * df_main$logFC_C1orf54_KD) +
  (weighted_fit_test_df_single_factor$logFC_lncRNA_KD * df_main$logFC_lncRNA_KD) #+
weighted_fit_test_df_single_factor$`weighted_fit_test_df_single_factor$residuals`


weighted_fit_test_real <-
  gam(y ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = weighted_fit_test_df_single_factor)

weighted_fit_test_real <-
  lm(y ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD), 
     data = weighted_fit_test_df_single_factor)
summary(weighted_fit_test_real)


#### make simple gam model use all genes
weighted_fit_test_df_single_factor <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD , 
      data = df_main_backup)
summary(weighted_fit_test_df_single_factor)

weighted_fit_test_df_single_factor <-
  lm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD, 
     data = raw_df)

df_new_model <-
  data.frame(beta_null = weighted_fit_test_df_single_factor$coefficients[1],
             VPS45_KD = weighted_fit_test_df_single_factor$coefficients[2] * 
               df_main_backup$logFC_VPS45_KD,
             C1orf54_KD = weighted_fit_test_df_single_factor$coefficients[4] * 
               df_main_backup$logFC_C1orf54_KD,
             lncRNA_KD = weighted_fit_test_df_single_factor$coefficients[3] * 
               df_main_backup$logFC_lncRNA_KD)


# df_new_model$residual <-
#   weighted_fit_test_df_single_factor$residuals

df_new_model$y <-
  df_new_model$beta_null +
  df_new_model$VPS45_KD +
  df_new_model$C1orf54_KD +
  df_new_model$lncRNA_KD 

df_new_model <-
  df_new_model[df_main_backup$Geneid %in% raw_df$Geneid, ]


df_new_model <-
  df_new_model[order(abs(df_new_model$y), decreasing = T), ]
# df_new_model <-
#   df_new_model[1:1267, ]

df_new_model <-
  df_new_model[1:1000, ]



test_results <-
  lm(y ~ 
       VPS45_KD +
       C1orf54_KD +
       lncRNA_KD +
       (VPS45_KD * C1orf54_KD) +
       (VPS45_KD * lncRNA_KD) +
       (C1orf54_KD * lncRNA_KD),
     data = df_new_model)
summary(test_results)

# save.image(file = "Regression_09Jan2023.RData")

(weighted_fit_test_df_single_factor$coefficients[2] * raw_df$logFC_VPS45_KD) +
  (weighted_fit_test_df_single_factor$coefficients[4] * raw_df$logFC_C1orf54_KD) +
  (weighted_fit_test_df_single_factor$coefficients[3] * raw_df$logFC_lncRNA_KD)
