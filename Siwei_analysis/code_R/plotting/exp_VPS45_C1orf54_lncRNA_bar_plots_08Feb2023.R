# Siwei 08 Feb 2023
# make expression bar plot of C1orf54, VPS45, lncRNA

# init
library(readxl)
library(ggplot2)
library(RColorBrewer)

df_to_plot <-
  read_excel("C1orf54_VPS45_lncRNA_bar_08Feb2023.xlsx")



ggplot(df_to_plot,
       aes(x = Days,
           y = Expression,
           colour = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone)),
              position = position_dodge2(width = -1)) +
  # position_dodge2(width = 2) +
  theme_classic()

##### plot by-line

df_to_plot <-
  read_excel("CD12_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "lncRNA")

ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                          name = "Dark2"),
                      name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 6,
                                          name = "Set2"),
                      name = "Clone") +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD12, lncRNA") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())

t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


## C1orf54
df_to_plot <-
  read_excel("CD12_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "C1orf54")

ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 6,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD12, C1orf54") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())

t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


## VPS45_all_transcripts
df_to_plot <-
  read_excel("CD12_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "VPS45_all_transcripts")

ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 6,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD12, VPS45_all") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())


t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


## VPS45_ENST369130
df_to_plot <-
  read_excel("CD12_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "VPS45_ENST369130")
ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 6,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD12, VPS45_369130") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())

t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


### CD_11

df_to_plot <-
  read_excel("CD11_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "lncRNA")

ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 7,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD11, lncRNA") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())

t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


## C1orf54
df_to_plot <-
  read_excel("CD11_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "C1orf54")

ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 7,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD11, C1orf54") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())

t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


## VPS45_all_transcripts
df_to_plot <-
  read_excel("CD11_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "VPS45_all_transcripts")

ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 7,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD11, VPS45_all") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())


t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)


## VPS45_ENST369130
df_to_plot <-
  read_excel("CD11_exp_relative_qPCR_table_4_plot.xlsx",
             sheet = "VPS45_ENST369130")
ggplot(df_to_plot,
       aes(x = Genotype,
           y = Rel_Exp_level,
           fill = as.factor(Genotype))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = as.factor(Clone_genotype)),
              width = 0.1,
              size = 2) +
  ylab("Relative expression level") +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2"),
                    name = "Genotype") +
  scale_colour_manual(values = brewer.pal(n = 7,
                                          name = "Set2"),
                      name = "Clone") +
  ylim(0, max(df_to_plot$Rel_Exp_level) * 1.2) +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("CD11, VPS45_369130") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank())

t.test(x = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "AA"],
       y = df_to_plot$Rel_Exp_level[df_to_plot$Genotype %in% "GG"],
       var.equal = F, alternative = "g", paired = F)

