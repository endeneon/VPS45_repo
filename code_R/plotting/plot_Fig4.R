# Siwei 10 Sept 2021
# plot Fig 4, S8-

# init
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(readxl)

##
# plot Fig4 PSD95_area_density

df_plot <- read_excel("Fig4_data.xlsx",
                      sheet = 1)

ggplot(df_plot,
       aes(x = Genotype,
           y = Value,
           fill = Genotype)) +
  geom_boxplot(outlier.alpha = 0,
               position = "dodge2",
               notch = F,
               size = 0.75) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               # position = "identity",
               # position = position_dodge(0.75),
               dotsize = 0.8,
               binwidth = max(df_plot$Value) / 50) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # ylim(0, 12) +
  xlab("") +
  ylab("PSD-95 puncta density (um-1)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))

##
# plot Fig4 SYN1_area_density

df_plot <- read_excel("Fig4_data.xlsx",
                      sheet = 2)

ggplot(df_plot,
       aes(x = Genotype,
           y = Value,
           fill = Genotype)) +
  geom_boxplot(outlier.alpha = 0,
               position = "dodge2",
               notch = F,
               size = 0.75) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               # position = "identity",
               # position = position_dodge(0.75),
               dotsize = 0.8,
               binwidth = max(df_plot$Value) / 50) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # ylim(0, 12) +
  xlab("") +
  ylab("SYN1 puncta density (um-2)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))


##
# plot Fig4 colocalisation_area_density

df_plot <- read_excel("Fig4_data.xlsx",
                      sheet = 3)

ggplot(df_plot,
       aes(x = Genotype,
           y = Value,
           fill = Genotype)) +
  geom_boxplot(outlier.alpha = 0,
               position = "dodge2",
               notch = F,
               size = 0.75) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               # position = "identity",
               # position = position_dodge(0.75),
               dotsize = 0.8,
               binwidth = max(df_plot$Value) / 50) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # ylim(0, 12) +
  xlab("") +
  ylab("Colocalisation puncta density (um-2)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))
