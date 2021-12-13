# Siwei 29 Oct 2020
# plot Fig S5

# init
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(readxl)

## import data
VPS45_qPCR <- read_excel("VPS45_qPCR.xlsx",
                         sheet = "VPS45")
lncRNA_qPCR <- read_excel("VPS45_qPCR.xlsx",
                          sheet = "lncRNA")
LINC00869_qPCR <- read_excel("VPS45_qPCR.xlsx",
                             sheet = "LINC00869")
C1orf54_qPCR <- read_excel("VPS45_qPCR.xlsx",
                           sheet = "C1orf54")

## t test
ttest_input <- read_excel("VPS45_qPCR.xlsx",
                         sheet = "VPS45")
ttest_input <- read_excel("VPS45_qPCR.xlsx",
                          sheet = "lncRNA")
ttest_input <- read_excel("VPS45_qPCR.xlsx",
                             sheet = "LINC00869")
ttest_input <- read_excel("VPS45_qPCR.xlsx",
                           sheet = "C1orf54")


t.test(x = ttest_input$`Relative Exp.`[ttest_input$Genotype == "AA"],
       y = ttest_input$`Relative Exp.`[ttest_input$Genotype == "GG"], 
       alternative = "t",
       paired = F,
       var.equal = F)


## plotting

ggplot(VPS45_qPCR,
       aes(x = Genotype,
           y = `Relative Exp.` * 1e4,
           fill = Genotype),
       shape = 0) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(stat = "identity",
              # position = "jitter",
              width = 0.04) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  ylim(0, 12) +
  ylab("VPS45 Relative Exp. x 1e-4") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(LINC00869_qPCR,
       aes(x = Genotype,
           y = `Relative Exp.` * 1e6,
           fill = Genotype),
       shape = 0) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(stat = "identity",
              # position = "jitter",
              width = 0.04) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  ylim(0, 12) +
  ylab("LINC00869 Relative Exp. x 1e-6") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(lncRNA_qPCR,
       aes(x = Genotype,
           y = `Relative Exp.` * 1e3,
           fill = Genotype),
       shape = 0) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(stat = "identity",
              # position = "jitter",
              width = 0.04) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  ylim(0, 3) +
  ylab("lncRNA C. Relative Exp. x 1e-3") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(C1orf54_qPCR,
       aes(x = Genotype,
           y = `Relative Exp.` * 1000,
           fill = Genotype),
       shape = 0) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(stat = "identity",
              # position = "jitter",
              width = 0.04) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # ylim(0, 12) +
  ylab("C1orf54 Relative Exp. x 1e-3") +
  theme_classic() +
  theme(legend.position = "none")

#### S5c
#### Use data from Fig 3b

ggplot(Fig3b_data, aes(x = `log2_FC_AA/GG`,
                       y = `log2_FC_AG/GG`)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = F) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ggtitle("Coefficient = 1.16,\nSpearman's Rho = 0.602, P = 0.00751") +
  theme_bw() +
  theme(axis.text = element_text(size = 14))

cor.test(x = Fig3b_data$`log2_FC_AG/GG`,
         y = Fig3b_data$`log2_FC_AA/GG`,
         method = "spearman")
lm.fit(
   x = as.matrix(Fig3b_data$`log2_FC_AG/GG`),
   y = as.matrix(Fig3b_data$`log2_FC_AA/GG`))
