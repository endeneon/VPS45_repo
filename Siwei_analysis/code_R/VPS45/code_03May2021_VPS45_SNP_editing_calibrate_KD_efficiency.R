# Siwei 03 May 2021
# Regression analysis on VPS45 A/G SNP editing vs DEG genes as Xin suggested
# use simple additive model, three explanatory variables (x)
# of lcRNA, vps45, and c1orf54
# try nonlinear variable model
# try to calibrate gene expression by scaling up/down gene KD efficiency
# note AA is the "baseline" in this dataset


# init

library(readr)
library(readxl)

library(MASS)
library(mgcv)

library(ppcor)
library(pscl)

library(ggplot2)
library(RColorBrewer)

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

raw_df <- raw_df[raw_df$FDR < 0.05, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_VPS45_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_lncRNA_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)




VPS45_GO_term_enriched_list <- read_excel("VPS45_GO_term_enriched_list.xlsx")
raw_df <- raw_df[raw_df$Geneid %in% VPS45_GO_term_enriched_list$Geneid, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = raw_df)

summary(weighted_fit_df)

# try nonlinear model
weighted_fit_df <-
  gam(logFC.groupGG ~ s(logFC_VPS45_KD) + s(logFC_lncRNA_KD) + s(logFC_C1orf54_KD),
      data = raw_df)

summary(weighted_fit_df)


# use heatmap FC enriched list
raw_df <- df_main

VPS45_heatmap_FC_enriched_list <- read_excel("VPS45_heatmap_FC_enriched_list.xlsx")
raw_df <- raw_df[raw_df$Gene_Symbol %in% VPS45_heatmap_FC_enriched_list$Gene_Symbol, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = raw_df)

summary(weighted_fit_df)

# try nonlinear (mixed) model
weighted_fit_df <-
  gam(logFC.groupGG ~ s(logFC_VPS45_KD) + s(logFC_lncRNA_KD) + s(logFC_C1orf54_KD),
      data = raw_df)

summary(weighted_fit_df)


##### test higher interactions #####

# use DEG gene list (N = 1267)
raw_df <- df_main

VPS45_GO_term_enriched_list <- read_excel("VPS45_GO_term_enriched_list.xlsx")
raw_df <- raw_df[raw_df$Geneid %in% VPS45_GO_term_enriched_list$Geneid, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ (logFC_VPS45_KD * logFC_lncRNA_KD) + logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_VPS45_KD + (logFC_lncRNA_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ (logFC_VPS45_KD * logFC_C1orf54_KD) + logFC_lncRNA_KD ,
      data = raw_df)
summary(weighted_fit_df)

# try nonlinear (mixed) model
weighted_fit_df <-
  gam(logFC.groupGG ~ s(logFC_VPS45_KD) + s(logFC_lncRNA_KD) + s(logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)


# use heatmap FC enriched list
raw_df <- df_main

VPS45_heatmap_FC_enriched_list <- read_excel("VPS45_heatmap_FC_enriched_list.xlsx")
raw_df <- raw_df[raw_df$Gene_Symbol %in% VPS45_heatmap_FC_enriched_list$Gene_Symbol, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

# try linear model
weighted_fit_df <-
  gam(logFC.groupGG ~ (logFC_VPS45_KD * logFC_lncRNA_KD) + logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_VPS45_KD + (logFC_lncRNA_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ (logFC_VPS45_KD * logFC_C1orf54_KD) + logFC_lncRNA_KD ,
      data = raw_df)
summary(weighted_fit_df)

# try nonlinear (mixed) model
weighted_fit_df <-
  gam(logFC.groupGG ~ s(logFC_VPS45_KD) + s(logFC_lncRNA_KD) + s(logFC_C1orf54_KD),
      data = raw_df)

summary(weighted_fit_df)


##### try high interaction linear model (complex)
# use all gene list
raw_df <- df_main

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  glm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +
        (logFC_lncRNA_KD * logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  lm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  lm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
     data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        # logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        s(logFC_lncRNA_KD * logFC_C1orf54_KD) +
        s(logFC_VPS45_KD * logFC_lncRNA_KD) +
        s(logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        (logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_lncRNA_KD * logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD * logFC_lncRNA_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD * logFC_C1orf54_KD,
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)


# use DEG gene list (N = 1267)
raw_df <- df_main

VPS45_GO_term_enriched_list <- read_excel("VPS45_GO_term_enriched_list.xlsx")
raw_df <- raw_df[raw_df$Geneid %in% VPS45_GO_term_enriched_list$Geneid, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  glm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)

# use heatmap FC enriched list
raw_df <- df_main
raw_df <- raw_df[raw_df$FDR < 0.05, ]

VPS45_heatmap_FC_enriched_list <- read_excel("VPS45_heatmap_FC_enriched_list.xlsx")
raw_df <- raw_df[raw_df$Gene_Symbol %in% VPS45_heatmap_FC_enriched_list$Gene_Symbol, ]
raw_df <- raw_df[!duplicated(raw_df$Geneid), ]

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)

weighted_fit_df <-
  glm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
      data = raw_df)
summary(weighted_fit_df)


weighted_fit_df <-
  gam(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD * logFC_C1orf54_KD),
      data = raw_df,
      gamma = 1)
summary(weighted_fit_df)
coef(weighted_fit_df)

weighted_fit_df <-
  glm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD) 
      data = raw_df)
summary(weighted_fit_df)
coef(weighted_fit_df)

save.image(file = "../R_synergy_analysis/VPS45_correlation_model_raw.RData")

load(file = "../R_synergy_analysis/VPS45_correlation_model_raw.RData")

weighted_fit_df <-
  gam(logFC.groupGG ~ 
        (logFC_lncRNA_KD * logFC_C1orf54_KD),
      data = raw_df)

weighted_fit_df <-
  lm(logFC.groupGG ~ 
        logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD), 
      data = raw_df)

weighted_fit_df <-
  lm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD), 
     data = raw_df)
summary(weighted_fit_df)

weighted_fit_df <-
  lm(logFC.groupGG ~ 
       # logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       # (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       # (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD), 
     data = raw_df)
summary(weighted_fit_df)

test_df <-
  data.frame(y = raw_df$logFC.groupGG,
             x = raw_df$logFC_VPS45_KD * raw_df$logFC_C1orf54_KD)
weighted_fit_df <-
  lm(y ~ x,
     data = test_df)
summary(weighted_fit_df)

weighted_fit_df$coefficients[4]

lm(logFC)

#### Try plots
library(ggplot2)

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_VPS45_KD)) +
  geom_point() +
  stat_smooth() +
  theme_classic()

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_VPS45_KD * logFC_C1orf54_KD)) +
  geom_point(size = 0.5) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  stat_smooth(method = "glm") +
  theme_bw()

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_VPS45_KD * logFC_lncRNA_KD)) +
  geom_point(size = 0.5) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  stat_smooth(method = "glm") +
  theme_bw()

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_C1orf54_KD * logFC_lncRNA_KD)) +
  geom_point(size = 0.5) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  stat_smooth(method = "glm") +
  theme_bw()

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_C1orf54_KD)) +
  geom_point(size = 0.5) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  stat_smooth(method = "glm") +
  theme_bw()

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_VPS45_KD)) +
  geom_point(size = 0.5) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  stat_smooth(method = "glm") +
  theme_bw()

ggplot(raw_df,
       aes(y = logFC.groupGG,
           x = logFC_lncRNA_KD)) +
  geom_point(size = 0.5) +
  xlim(-2, 2) +
  ylim(-2, 2) +
  stat_smooth(method = "glm") +
  theme_bw()

cor.test(x = raw_df$logFC.groupGG,
         y = raw_df$logFC_VPS45_KD * raw_df$logFC_C1orf54_KD,
         method = "pearson",
         alternative = "t")

cor.test(x = raw_df$logFC.groupGG,
         y = raw_df$logFC_VPS45_KD * raw_df$logFC_lncRNA_KD,
         method = "spearman",
         alternative = "t")

cor.test(x = df_main$logFC.groupGG,
         y = df_main$logFC_VPS45_KD * df_main$logFC_C1orf54_KD,
         method = "pearson", 
         alternative = "t")

cor.test(x = raw_df$logFC.groupGG,
         y = raw_df$logFC_VPS45_KD * raw_df$logFC_C1orf54_KD,
         method = "pearson",
         alternative = "t")

# C1orf54
t.test(x = c(1.3905383380953873, 1.4150802887640894),
       y = c(1.1364051189626745, 1.099328239300169), 
       alternative = "t", paired = F, var.equal = F)

# VPS45
t.test(x = c(7.614072534907256, 7.543755098042264),
       y = c(7.107189988932039, 7.188836053027453), 
       alternative = "t", paired = F, var.equal = F)


### check the correlation of ground vs 3xKD
results_QLM_3xKD_vs_AA_sub <-
  results_QLM_3xKD_vs_AA
results_QLM_3xKD_vs_AA_sub$Geneid <-
  rownames(results_QLM_3xKD_vs_AA_sub)
results_QLM_3xKD_vs_AA_sub <-
  results_QLM_3xKD_vs_AA_sub[results_QLM_3xKD_vs_AA_sub$Geneid %in% df_main$Geneid, ]
results_QLM_3xKD_vs_AA_sub <-
  results_QLM_3xKD_vs_AA_sub[, c("Geneid", "logFC", "FDR")]
colnames(results_QLM_3xKD_vs_AA_sub)[2] <- "logFC_KDx3"

raw_df <-
  merge(x = df_main,
        y = results_QLM_3xKD_vs_AA_sub,
        by = "Geneid")

weighted_fit_df <-
  lm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD),
     data = raw_df)


weighted_fit_df <-
  glm(logFC.groupGG ~ 
       logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
       (logFC_lncRNA_KD * logFC_C1orf54_KD) +
       (logFC_VPS45_KD * logFC_lncRNA_KD) +
       (logFC_VPS45_KD * logFC_C1orf54_KD),
     data = raw_df)
# weighted_fit_df <-
#   gam(logFC.groupGG ~ 
#         # logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
#         # (logFC_lncRNA_KD * logFC_C1orf54_KD) +
#         # (logFC_VPS45_KD * logFC_lncRNA_KD) +
#         # (logFC_VPS45_KD * logFC_C1orf54_KD) +
#         (logFC_VPS45_KD + logFC_lncRNA_KD * logFC_C1orf54_KD),
      # data = raw_df)
summary(weighted_fit_df)
pR2(weighted_fit_df)['McFadden']
weighted_fit_df$coefficients
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      -0.018928   0.012842  -1.474  0.14076    
# logFC_VPS45_KD                   -0.140267   0.036360  -3.858  0.00012 ***
#   logFC_lncRNA_KD                  -0.148252   0.027311  -5.428 6.82e-08 ***
#   logFC_C1orf54_KD                 -0.009173   0.034812  -0.264  0.79220    
# logFC_lncRNA_KD:logFC_C1orf54_KD -0.320915   0.055478  -5.785 9.16e-09 ***
#   logFC_VPS45_KD:logFC_lncRNA_KD    0.164483   0.041114   4.001 6.68e-05 ***
#   logFC_VPS45_KD:logFC_C1orf54_KD   0.263763   0.060556   4.356 1.43e-05 ***

####
weighted_fit_df <-
  glm(logFC.groupGG ~ 
        # logFC_VPS45_KD #+ 
        # logFC_lncRNA_KD #+
        # logFC_C1orf54_KD #+
        # (logFC_lncRNA_KD * logFC_C1orf54_KD) #+
        # (logFC_VPS45_KD * logFC_lncRNA_KD) #+
        (logFC_VPS45_KD * logFC_C1orf54_KD)
      ,
      data = raw_df)
pR2(weighted_fit_df)['McFadden']






library(ppcor)

df_partial_test <-
  data.frame(cbind(raw_df$logFC.groupGG,
                   raw_df$logFC_VPS45_KD,
                   raw_df$logFC_lncRNA_KD,
                   raw_df$logFC_C1orf54_KD),
             stringsAsFactors = F)
colnames(df_partial_test) <-
  c("groupGG", "VPS45_KD", "lncRNA_KD", "C1orf54_KD")
spcor(df_partial_test,
      method = "pearson")

#### 3xKD vs AA needs to use its own GG vs AA
results_QLM_3x_all <-
  cbind(results_QLM_GG_vs_AA[, c("logFC", "PValue")],
        results_QLM_3xKD_vs_AA[, c("logFC", "PValue")])
colnames(results_QLM_3x_all) <-
  c("logFC_GGvsAA", "Pval_GGvsAA", 
    "logFC_3xKDvsAA", "Pval_3xGGvsAA")
results_QLM_3x_all$Geneid <-
  rownames(results_QLM_3xKD_vs_AA)
results_QLM_3x_all <-
  results_QLM_3x_all[results_QLM_3x_all$Geneid %in% raw_df$Geneid, ]

weighted_fit_df_3xKD <-
  lm(logFC_GGvsAA ~ 
       logFC_3xKDvsAA, 
     data = results_QLM_3x_all)
cor.test(x = results_QLM_3x_all$logFC_GGvsAA,
         y = results_QLM_3x_all$logFC_3xKDvsAA,
         method = "pearson")

cor.test(x = raw_df$logFC.groupGG,
         y = raw_df$logFC_VPS45_KD,
         method = "pearson", alternative = "t")

summary(weighted_fit_df_3xKD)
