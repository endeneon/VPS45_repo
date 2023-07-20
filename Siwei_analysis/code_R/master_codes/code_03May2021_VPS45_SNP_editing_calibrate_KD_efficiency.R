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
        # logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD +
        (logFC_lncRNA_KD * logFC_C1orf54_KD) +
        (logFC_VPS45_KD * logFC_lncRNA_KD) +
        (logFC_VPS45_KD * logFC_C1orf54_KD),
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
         method = "spearman")

# C1orf54
t.test(x = c(1.3905383380953873, 1.4150802887640894),
       y = c(1.1364051189626745, 1.099328239300169), 
       alternative = "t", paired = F, var.equal = F)

# VPS45
t.test(x = c(7.614072534907256, 7.543755098042264),
       y = c(7.107189988932039, 7.188836053027453), 
       alternative = "t", paired = F, var.equal = F)
