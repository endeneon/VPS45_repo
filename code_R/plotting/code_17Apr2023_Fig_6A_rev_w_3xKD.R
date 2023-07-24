# Siwei 14 Apr 2023
# Siwei 31 Aug 2021
# Re-plot Fig 6A w/ data of 3xKD

# init
library(ggplot2)
library(RColorBrewer)

library(readxl)
library(stringr)

# load data

## load VPS45_iso
df_to_plot <-
  read_excel("Fig6A_data_rev_17Apr2023.xlsx",
             sheet = 1)
sheet_title <-
  excel_sheets(path = "Fig6A_data_rev_17Apr2023.xlsx")[1]

## load lncRNA
df_to_plot <-
  read_excel("Fig6A_data_rev_17Apr2023.xlsx",
             sheet = 2)
sheet_title <-
  excel_sheets(path = "Fig6A_data_rev_17Apr2023.xlsx")[2]

## load C1orf54
df_to_plot <-
  read_excel("Fig6A_data_rev_17Apr2023.xlsx",
             sheet = 3)
sheet_title <-
  excel_sheets(path = "Fig6A_data_rev_17Apr2023.xlsx")[3]

## load VPS45_all
df_to_plot <-
  read_excel("Fig6A_data_rev_17Apr2023.xlsx",
             sheet = 4)
sheet_title <-
  excel_sheets(path = "Fig6A_data_rev_17Apr2023.xlsx")[4]

df_to_plot$`Genotype and Knockdown gene` <-
  factor(df_to_plot$`Genotype and Knockdown gene`,
         levels = c("AA_EGFP_KD",
                    "AA_VPS45_KD",
                    "AA_lncRNA_KD",
                    "AA_C1orf54_KD",
                    "AA_3x_KD",
                    "GG_EGFP_KD"))
# make the plot
ggplot(df_to_plot,
       aes(x = `Genotype and Knockdown gene`,
           y = `Relative Expression Value`,
           fill = `Genotype and Knockdown gene`)) +
  geom_boxplot(outlier.colour = NULL,
               outlier.shape = NA,
               outlier.size = 0,
               notch = F,
               width = 0.8) +
  geom_point() +
  stat_summary(fun = mean,
               geom ="point",
               shape = 4,
               size= 2) +
  xlab(sheet_title) +
  scale_fill_discrete(name = "Genotypes and \nKnockdown Genes",
                      type = c("#1B9E77",
                               "#E7298A",
                               "#D95F02",
                               "#666666",
                               "yellow3",
                               "#7570B3")) +
  ylim(0, max(df_to_plot$`Relative Expression Value`)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   vjust = 0,
                                   hjust = 0))


# t = 11.362, df = 5.4689, p-value = 5.202e-05
t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP_KD"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_3x_KD"],
       var.equal = F, alternative = "t", paired = F)

## load lncRNA
df_to_plot <-
  read_excel("Fig6A_data.xlsx",
             sheet = 2)
# remove all "_KD"
df_to_plot$`Genotype and Knockdown gene` <-
  str_replace_all(string = df_to_plot$`Genotype and Knockdown gene`,
                  pattern = "_KD$",
                  replacement = "")


df_to_plot$`Genotype and Knockdown gene` <-
  factor(df_to_plot$`Genotype and Knockdown gene`,
         levels = c("AA_EGFP",
                    "AA_VPS45",
                    "AA_lncRNA",
                    "AA_C1orf54",
                    "GG_EGFP"))
# make the plot
ggplot(df_to_plot,
       aes(x = `Genotype and Knockdown gene`,
           y = `Relative Expression Value`,
           fill = `Genotype and Knockdown gene`)) +
  geom_boxplot(outlier.colour = NULL,
               outlier.shape = NA,
               outlier.size = 0,
               notch = F,
               width = 0.8) +
  stat_summary(fun = mean,
               geom ="point",
               shape = 4,
               size= 2) +
  xlab("AC244033.2") +
  scale_fill_discrete(name = "Genotypes and \nKnockdown Genes",
                      type = c("#1B9E77",
                               "#E7298A",
                               "#D95F02",
                               "#666666",
                               "#7570B3")) +
  ylim(0, max(df_to_plot$`Relative Expression Value`)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   vjust = 0,
                                   hjust = 0))
## run t.test
t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_VPS45"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_lncRNA"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_C1orf54"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "GG_EGFP"],
       var.equal = F, alternative = "t", paired = F)

 ## load C1orf54
df_to_plot <-
  read_excel("Fig6A_data.xlsx",
             sheet = 3)
# remove all "_KD"
df_to_plot$`Genotype and Knockdown gene` <-
  str_replace_all(string = df_to_plot$`Genotype and Knockdown gene`,
                  pattern = "_KD$",
                  replacement = "")


df_to_plot$`Genotype and Knockdown gene` <-
  factor(df_to_plot$`Genotype and Knockdown gene`,
         levels = c("AA_EGFP",
                    "AA_VPS45",
                    "AA_lncRNA",
                    "AA_C1orf54",
                    "GG_EGFP"))
# make the plot
ggplot(df_to_plot,
       aes(x = `Genotype and Knockdown gene`,
           y = `Relative Expression Value`,
           fill = `Genotype and Knockdown gene`)) +
  geom_boxplot(outlier.colour = NULL,
               outlier.shape = NA,
               outlier.size = 0,
               notch = F,
               width = 0.8) +
  stat_summary(fun = mean,
               geom ="point",
               shape = 4,
               size= 2) +
  xlab("C1ORF54") +
  scale_fill_discrete(name = "Genotypes and \nKnockdown Genes",
                      type = c("#1B9E77",
                               "#E7298A",
                               "#D95F02",
                               "#666666",
                               "#7570B3")) +
  ylim(0, max(df_to_plot$`Relative Expression Value`)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   vjust = 0,
                                   hjust = 0))

## run t.test
t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_VPS45"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_lncRNA"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_C1orf54"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "GG_EGFP"],
       var.equal = F, alternative = "t", paired = F)

 ## load VPS45_all
df_to_plot <-
  read_excel("Fig6A_data.xlsx",
             sheet = 4)
# remove all "_KD"
df_to_plot$`Genotype and Knockdown gene` <-
  str_replace_all(string = df_to_plot$`Genotype and Knockdown gene`,
                  pattern = "_KD$",
                  replacement = "")


df_to_plot$`Genotype and Knockdown gene` <-
  factor(df_to_plot$`Genotype and Knockdown gene`,
         levels = c("AA_EGFP",
                    "AA_VPS45",
                    "AA_lncRNA",
                    "AA_C1orf54",
                    "GG_EGFP"))
# make the plot
ggplot(df_to_plot,
       aes(x = `Genotype and Knockdown gene`,
           y = `Relative Expression Value`,
           fill = `Genotype and Knockdown gene`)) +
  geom_boxplot(outlier.colour = NULL,
               outlier.shape = NA,
               outlier.size = 0,
               notch = F,
               width = 0.8) +
  stat_summary(fun = mean,
               geom ="point",
               shape = 4,
               size= 2) +
  xlab("VPS45 all isoforms") +
  scale_fill_discrete(name = "Genotypes and \nKnockdown Genes",
                      type = c("#1B9E77",
                               "#E7298A",
                               "#D95F02",
                               "#666666",
                               "#7570B3")) +
  ylim(0, max(df_to_plot$`Relative Expression Value`)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   vjust = 0,
                                   hjust = 0))

## run t.test
t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_VPS45"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_lncRNA"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_C1orf54"],
       var.equal = F, alternative = "t", paired = F)

t.test(x = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "AA_EGFP"],
       y = df_to_plot$`Relative Expression Value`[df_to_plot$`Genotype and Knockdown gene` %in% "GG_EGFP"],
       var.equal = F, alternative = "t", paired = F)
