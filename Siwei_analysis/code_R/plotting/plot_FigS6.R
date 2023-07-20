# Siwei 10 Jan 2021
# plot Fig S6-S7

# init
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(readxl)

# plot FigS7A
S7A_data <- read_excel("VPS45_bulk_vs_CROPSeq_FDR0.2.xlsx",
                       sheet = "full_list")
# filter for FDR < 0.2 in both sets
S7A_data <- S7A_data[(S7A_data$FDR_CRISPR < 0.2) & (S7A_data$FDR_CROP < 0.2), ]

ggplot(S7A_data,
       aes(x = logFC.groupGG,
           y = FC_CROPSeq)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlim(-1, 1) +
  ylim(-2, 2) +
  xlab("Log2 fold change in bulk RNA-seq (SNP editing)") +
  ylab("Log2 fold change in VPS45 CROP-seq") +
  stat_smooth(method = "lm",
              se = F) +
  ggtitle("Spearman's rho = 0.314, P = 7.724x10-5") +
  theme_classic()

cor.test(S7A_data$logFC.groupGG,
         S7A_data$FC_CROPSeq,
         method = "s",
         alternative = "t")
