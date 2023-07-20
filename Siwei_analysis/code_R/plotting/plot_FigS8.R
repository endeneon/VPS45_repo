# Siwei 10 Sept 2021
# plot Fig 4, S8-

# init
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(readxl)

##
# plot FigS8 spine_protrusion

df_plot <- read_excel("FigS8_data.xlsx",
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
  ylab("Protrusions (um-1)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))

##
# plot FigS8 PSD95 puncta area

df_plot <- read_excel("FigS8_data.xlsx",
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
  ylab("PSD-95 puncta area (um2)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))


##
# plot FigS8 SYN1 puncta area

df_plot <- read_excel("FigS8_data.xlsx",
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
  ylab("SYN1 puncta area (um2)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))

##
# plot FigS8 colocalisation puncta area

df_plot <- read_excel("FigS8_data.xlsx",
                      sheet = 4)

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
  ylab("Colocalisation puncta area (um2)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))
