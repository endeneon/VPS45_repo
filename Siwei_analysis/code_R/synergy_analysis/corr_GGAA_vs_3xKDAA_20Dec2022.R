# Siwei 20 Dec 2022
# calc the correlation of GG/AA vs 3xKD/AA
# may need to pick the top x FDR sigif genes

# init
library(MASS)
library(mgcv)

library(ggplot2)
library(RColorBrewer)



# load data

# load("~/NVME/VPS45/R_synergy_analysis/df_results_GG_3xKD_AA.RData")
load("~/NVME/VPS45/R_synergy_analysis/df_results_GG_3xKD_AA_no_sg2.RData")
load("~/NVME/VPS45/R_synergy_analysis/VPS45_correlation_model_raw.RData")

results_QLM_3xKD_vs_AA$Geneid <- rownames(results_QLM_3xKD_vs_AA)
results_QLM_3xKD_vs_AA <-
  results_QLM_3xKD_vs_AA[, c("Geneid", "logFC", "FDR")]
results_QLM_GG_vs_AA$Geneid <- rownames(results_QLM_GG_vs_AA)
results_QLM_GG_vs_AA <-
  results_QLM_GG_vs_AA[, c("Geneid", "logFC", "FDR")]

raw_df <-
  merge(x = raw_df,
        y = results_QLM_3xKD_vs_AA,
        by = "Geneid")
colnames(raw_df)[c(ncol(raw_df) - 1,
                   ncol(raw_df))] <-
  c("logFC_3xKD_AA", "FDR_3xKD_AA")

raw_df <-
  merge(x = raw_df,
        y = results_QLM_GG_vs_AA,
        by = "Geneid")
colnames(raw_df)[c(ncol(raw_df) - 1,
                   ncol(raw_df))] <-
  c("logFC_3xGG_vs_AA", "FDR_3xGG_AA")

weighted_fit_df <-
  gam(logFC.groupGG ~ 0 +
        # logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +
        logFC_3xGG_vs_AA +
        logFC_3xKD_AA,
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xGG_vs_AA ~ 0 +
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 0 +
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)



weighted_fit_df <-
  gam(logFC.groupGG ~ 0 +
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
        # (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD) +
        # logFC_3xGG_vs_AA +
        # logFC_3xKD_AA,
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xGG_vs_AA ~ 0 +
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        logFC_3xKD_AA,
        # (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        # (logFC_VPS45_KD * logFC_lncRNA_KD) +
        # (logFC_VPS45_KD * logFC_C1orf54_KD) +
        # (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xGG_vs_AA ~ 0 +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +
        # (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD) +
        logFC_3xKD_AA,
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xGG_vs_AA ~ 0 +
        logFC_3xKD_AA,
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xKD_AA ~ 0 +
        # (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        # (logFC_VPS45_KD * logFC_lncRNA_KD) +
        # (logFC_VPS45_KD * logFC_C1orf54_KD) +
        (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xKD_AA ~ 0 +
        # (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        # (logFC_VPS45_KD * logFC_lncRNA_KD) +
        # (logFC_VPS45_KD * logFC_C1orf54_KD) +
        (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xKD_AA ~ 0 +
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,# +
        # (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        # (logFC_VPS45_KD * logFC_lncRNA_KD) +
        # (logFC_VPS45_KD * logFC_C1orf54_KD) +# +
        # (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  gam(logFC_3xKD_AA ~ 0 +
        # logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,# +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +# +
        (logFC_lncRNA_KD* logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)



cor.test(x = results_QLM_GG_vs_AA$logFC,
         y = results_QLM_3xKD_vs_AA$logFC,
         alternative = "t",
         method = "s")
