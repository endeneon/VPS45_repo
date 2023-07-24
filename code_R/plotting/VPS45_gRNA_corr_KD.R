# Siwei 01 Feb 2023
# plot VPS45 gRNA1/2/3 CROP-seq vs KD results

# init
library(readxl)
library(ggplot2)

VPS45_correlation_spearman_calc <-
  read_excel("VPS45_correlation_spearman_calc.xlsx",
             sheet = 1)

df_to_plot <-
  data.frame(gRNA = "gRNA",
             gene_KD = "gene_KD",
             rho = "1",
             p.value = "1")

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...2,
           y = VPS45_correlation_spearman_calc$VPS45_KD,
           method = "spearman",
           alternative = "t",
           na.action = "na.omit")
df_to_plot[1, ] <-
  c("VPS45_gRNA_1", "VPS45_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)


df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...5,
         y = VPS45_correlation_spearman_calc$lncRNA_KD,
         method = "spearman",
         alternative = "t",
         na.action = "na.omit")
df_to_plot[2, ] <-
  c("VPS45_gRNA_1", "lncRNA_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...8,
         y = VPS45_correlation_spearman_calc$C1orf54_KD,
         method = "spearman",
         alternative = "t",
         exact = F,
         na.action = "na.omit")
df_to_plot[3, ] <-
  c("VPS45_gRNA_1", "C1orf54_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$logFC,
         y = VPS45_correlation_spearman_calc$Triple_KD,
         method = "spearman",
         alternative = "t",
         # exact = F,
         na.action = "na.omit")
df_to_plot[4, ] <-
  c("VPS45_gRNA_1", "KDx3",
    df_cor$estimate,
    p.value = df_cor$p.value)

### gRNA 2
VPS45_correlation_spearman_calc <-
  read_excel("VPS45_correlation_spearman_calc.xlsx",
             sheet = 2)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...2,
           y = VPS45_correlation_spearman_calc$VPS45_KD,
           method = "spearman",
           alternative = "t",
           exact = F,
           na.action = "na.omit")
df_to_plot[5, ] <-
  c("VPS45_gRNA_2", "VPS45_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)


df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...5,
           y = VPS45_correlation_spearman_calc$lncRNA_KD,
           method = "spearman",
           alternative = "t",
           na.action = "na.omit")
df_to_plot[6, ] <-
  c("VPS45_gRNA_2", "lncRNA_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...8,
           y = VPS45_correlation_spearman_calc$C1orf54_KD,
           method = "spearman",
           alternative = "t",
           exact = F,
           na.action = "na.omit")
df_to_plot[7, ] <-
  c("VPS45_gRNA_2", "C1orf54_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$logFC,
           y = VPS45_correlation_spearman_calc$Triple_KD,
           method = "spearman",
           alternative = "t",
           # exact = F,
           na.action = "na.omit")
df_to_plot[8, ] <-
  c("VPS45_gRNA_2", "KDx3",
    df_cor$estimate,
    p.value = df_cor$p.value)


### gRNA 3
VPS45_correlation_spearman_calc <-
  read_excel("VPS45_correlation_spearman_calc.xlsx",
             sheet = 3)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...2,
           y = VPS45_correlation_spearman_calc$VPS45_KD,
           method = "spearman",
           alternative = "t",
           # exact = F,
           na.action = "na.omit")
df_to_plot[9, ] <-
  c("VPS45_gRNA_3", "VPS45_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)


df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...5,
           y = VPS45_correlation_spearman_calc$lncRNA_KD,
           method = "spearman",
           alternative = "t",
           na.action = "na.omit")
df_to_plot[10, ] <-
  c("VPS45_gRNA_3", "lncRNA_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_gRNA1...8,
           y = VPS45_correlation_spearman_calc$C1orf54_KD,
           method = "spearman",
           alternative = "t",
           # exact = F,
           na.action = "na.omit")
df_to_plot[11, ] <-
  c("VPS45_gRNA_3", "C1orf54_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$logFC,
           y = VPS45_correlation_spearman_calc$Triple_KD,
           method = "spearman",
           alternative = "t",
           # exact = F,
           na.action = "na.omit")
df_to_plot[12, ] <-
  c("VPS45_gRNA_3", "KDx3",
    df_cor$estimate,
    p.value = df_cor$p.value)

### CRISPR
VPS45_correlation_spearman_calc <-
  read_excel("VPS45_correlation_spearman_calc.xlsx",
             sheet = 4)

df_cor <-
  cor.test(x = 0 - VPS45_correlation_spearman_calc$VPS45_CRISPR...2,
           y = VPS45_correlation_spearman_calc$VPS45_KD,
           method = "spearman",
           alternative = "t",
           # exact = F,
           na.action = "na.omit")
df_to_plot[13, ] <-
  c("VPS45_CRISPR_GtoA", "VPS45_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)


df_cor <-
  cor.test(x = 0 - VPS45_correlation_spearman_calc$VPS45_CRISPR...5,
           y = VPS45_correlation_spearman_calc$lncRNA_KD,
           method = "spearman",
           alternative = "t",
           exact = F,
           na.action = "na.omit")
df_to_plot[14, ] <-
  c("VPS45_CRISPR_GtoA", "lncRNA_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = 0 - VPS45_correlation_spearman_calc$VPS45_CRISPR...8,
           y = VPS45_correlation_spearman_calc$C1orf54_KD,
           method = "spearman",
           alternative = "t",
           exact = F,
           na.action = "na.omit")
df_to_plot[15, ] <-
  c("VPS45_CRISPR_GtoA", "C1orf54_KD",
    df_cor$estimate,
    p.value = df_cor$p.value)

df_cor <-
  cor.test(x = VPS45_correlation_spearman_calc$VPS45_CRISPR...11,
           y = VPS45_correlation_spearman_calc$Triple_KD,
           method = "spearman",
           alternative = "t",
           exact = F,
           na.action = "na.omit")
df_to_plot[16, ] <-
  c("VPS45_CRISPR_GtoA", "KDx3",
    df_cor$estimate,
    p.value = df_cor$p.value)

ggplot(df_to_plot,
       aes(x = factor(gRNA,
                      levels = c("VPS45_gRNA_1",
                                 "VPS45_gRNA_2",
                                 "VPS45_gRNA_3",
                                 "VPS45_CRISPR_GtoA")),
           y = factor(gene_KD,
                      levels = c("VPS45_KD",
                                 "lncRNA_KD",
                                 "C1orf54_KD",
                                 "KDx3")),
           fill = as.numeric(rho),
           size = 0 - log10(as.numeric(p.value) + 5e-14))) +
  geom_point(shape = 21) +
  scale_radius(name = "-log10P",
               range = c(1, 10),
               limits = c(0,
                          max(0 - log10(as.numeric(df_to_plot$p.value) + 5e-14)))) +
  scale_fill_gradient2(low = "#4444CC",
                       mid = "white",
                       high = "#CC4444",
                       midpoint = 0,
                       # limits = c(-0.5, 0.5),
                       name = "Spearman's Rho") +
  scale_y_discrete(limits = rev) +
  xlab("") +
  ylab("shRNA knockdown") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   size = 12),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), "pt"))


save.image("VPS45_corr_02Feb2023.RData")

margin(1,1,1,1, unit = "pt")

df_to_plot <-
  df_to_plot[!(df_to_plot$gRNA %in% "VPS45_CRISPR_GtoA"), ]

CRISPR_AA_GG_vs_KDx3[CRISPR_AA_GG_vs_KDx3$Geneid %in% "ENSG00000369130", ]

sum(CRISPR_AA_GG_vs_KDx3$PValue...4 < 0.2)
sum(CRISPR_AA_GG_vs_KDx3$PValue...9 < 0.2)
sum(sum(CRISPR_AA_GG_vs_KDx3$PValue...4 < 0.2))