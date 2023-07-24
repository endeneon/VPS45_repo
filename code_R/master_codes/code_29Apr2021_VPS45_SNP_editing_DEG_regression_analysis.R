# Siwei 29 Apr 2021
# Regression analysis on VPS45 A/G SNP editing vs DEG genes as Xin suggested
# use simple additive model, three explanatory variables (x)
# of lcRNA, vps45, and c1orf54
# note AA is the "baseline" in this dataset


# init

library(readr)
library(readxl)

library(MASS)
library(mgcv)

library(ggplot2)
library(RColorBrewer)

# load data
raw_df <- read_excel("VPS45_SNP_editing_DEG_KD_regression.xlsx", 
                     sheet = "DGE_Overalp_w_KD")

raw_df$logFDR.vps45 <- 0 - log10(raw_df$FDR)


# try use VPS45 KD data series first, use rlm()

weighted_fit_df <- 
  rlm(logFC.groupGG ~ KD_vps45,
      data = raw_df,
      weights = raw_df$logFDR.vps45/sum(raw_df$logFDR.vps45))

weighted_fit_df <- 
  rlm(logFC.groupGG ~ KD_lncRNA,
      data = raw_df,
      weights = raw_df$logFDR.vps45/sum(raw_df$logFDR.vps45))

weighted_fit_df <- 
  rlm(logFC.groupGG ~ KD_c1orf54,
      data = raw_df,
      weights = raw_df$logFDR.vps45/sum(raw_df$logFDR.vps45))

# weighted_fit_df <- 
#   rlm(logFC.groupGG ~ KD_vps45,
#       data = raw_df)

qqnorm(weighted_fit_df$residuals,
       ylim = c(-2, 2))

####

## VPS45
# raw_df_005 <- raw_df[raw_df$FDR.vps45 < 0.05, ]

ggplot(raw_df, aes(x = KD_vps45, 
                       y = logFC.groupGG)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(colour = brewer.pal(4, "Dark2")[1],
             alpha = 0.5, 
             size = 0.5) +
  geom_smooth(method = "lm",
              colour = "magenta2") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  ggtitle(label = paste("rs2027349 A>G ~ VPS45 KD, beta =",
                        format(as.numeric(lm(logFC.groupGG ~ KD_vps45,
                                             data = raw_df)$coefficients[2]),
                               digits = 3),
                        sep = " ")) +
  theme_bw()

## lncRNA
ggplot(raw_df, aes(x = KD_lncRNA, 
                   y = logFC.groupGG)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color = brewer.pal(4, "Dark2")[2],
             alpha = 0.5, 
             size = 0.5) +
  geom_smooth(method = "lm",
              colour = "magenta2") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  ggtitle(label = paste("rs2027349 A>G ~ lncRNA KD, beta =",
                        format(as.numeric(lm(logFC.groupGG ~ KD_lncRNA,
                                             data = raw_df)$coefficients[2]),
                               digits = 3),
                        sep = " ")) +
  theme_bw()

## c1orf54
ggplot(raw_df, aes(x = KD_c1orf54, 
                   y = logFC.groupGG)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color = brewer.pal(4, "Dark2")[3],
             alpha = 0.5, 
             size = 0.5) +
  geom_smooth(method = "lm",
              colour = "magenta2") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  ggtitle(label = paste("rs2027349 A>G ~ c1orf54 KD, beta =",
                        format(as.numeric(lm(logFC.groupGG ~ KD_c1orf54,
                                             data = raw_df)$coefficients[2]),
                               digits = 3),
                        sep = " ")) +
  theme_bw()


lm(logFC.groupGG ~ KD_vps45 + KD_lncRNA + KD_c1orf54,
   data = raw_df)

### 
## VPS45
raw_df_005 <- raw_df[raw_df$FDR.vps45 < 0.05, ]

ggplot(raw_df_005, aes(x = KD_vps45, 
                   y = logFC.groupGG)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(colour = brewer.pal(4, "Dark2")[1],
             alpha = 0.5, 
             size = 0.5) +
  geom_smooth(method = "lm",
              colour = "magenta2") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  ggtitle(label = paste("DEG FDR< 0.05\n", 
                        "rs2027349 A>G ~ VPS45 KD, beta =",
                        format(as.numeric(lm(logFC.groupGG ~ KD_vps45,
                                             data = raw_df_005)$coefficients[2]),
                               digits = 3),
                        sep = " ")) +
  theme_bw()

## lncRNA
raw_df_005 <- raw_df[raw_df$FDR.lcRNA < 0.05, ]

ggplot(raw_df_005, aes(x = KD_lncRNA, 
                   y = logFC.groupGG)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color = brewer.pal(4, "Dark2")[2],
             alpha = 0.5, 
             size = 0.5) +
  geom_smooth(method = "lm",
              colour = "magenta2") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  ggtitle(label = paste("DEG FDR< 0.05\n", 
                        "rs2027349 A>G ~ lncRNA KD, beta =",
                        format(as.numeric(lm(logFC.groupGG ~ KD_lncRNA,
                                             data = raw_df_005)$coefficients[2]),
                               digits = 3),
                        sep = " ")) +
  theme_bw()

## c1orf54
raw_df_005 <- raw_df[raw_df$FDR.c1orf54 < 0.05, ]

ggplot(raw_df_005, aes(x = KD_c1orf54, 
                   y = logFC.groupGG)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color = brewer.pal(4, "Dark2")[3],
             alpha = 0.5, 
             size = 0.5) +
  geom_smooth(method = "lm",
              colour = "magenta2") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  ggtitle(label = paste("DEG FDR< 0.05\n", 
                        "rs2027349 A>G ~ c1orf54 KD, beta =",
                        format(as.numeric(lm(logFC.groupGG ~ KD_c1orf54,
                                             data = raw_df_005)$coefficients[2]),
                               digits = 3),
                        sep = " ")) +
  theme_bw()

weighted_fit_df <-
  lm(logFC.groupGG ~ KD_vps45 + KD_lncRNA + KD_c1orf54,
     data = raw_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ KD_vps45 + KD_lncRNA + KD_c1orf54,
      data = raw_df)

summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ s(KD_vps45) + s(KD_lncRNA) + s(KD_c1orf54),
      data = raw_df)

summary(weighted_fit_df)

cor.test(raw_df$logFC.groupGG,
         raw_df$KD_vps45)
cor.test(raw_df$logFC.groupGG,
         raw_df$KD_lncRNA)

qqnorm(weighted_fit_df$residuals)


##### use all genes #####

df_main <- read_excel("VPS45_AA_AG_GG_all_gene_table.xlsx", 
                      sheet = "VPS45_AA_AG_GG_all_gene_table_w")

df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,6)]
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_VPS45_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_VPS45_KD")
rm(df_to_merge)

df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_lncRNAKD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,6)]
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_lncRNA_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_lncRNA_KD")
rm(df_to_merge)

df_to_merge <- read_delim("Rapid_neuron_gene_exp_table_C1orf54KD_vs_EGFP_31Mar2021.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
df_to_merge <- df_to_merge[, c(1,2,6)]
colnames(df_to_merge) <- paste0(colnames(df_to_merge), "_C1orf54_KD")
df_main <- merge(df_main, df_to_merge,
                 by.x = "Geneid",
                 by.y = "Geneid_C1orf54_KD")
rm(df_to_merge)

##
raw_df <- df_main

weighted_fit_df <-
  gam(logFC.groupGG ~ logFC_VPS45_KD + logFC_lncRNA_KD + logFC_C1orf54_KD,
      data = raw_df)

summary(weighted_fit_df)

weighted_fit_df <-
  gam(logFC.groupGG ~ s(logFC_VPS45_KD) + s(logFC_lncRNA_KD) + s(logFC_C1orf54_KD),
      data = raw_df)

summary(weighted_fit_df)
