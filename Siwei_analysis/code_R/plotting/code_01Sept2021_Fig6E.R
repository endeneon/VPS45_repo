# Siwei 01 Sept 2021
# Re-plot Fig 6E

# init
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(colorspace)

library(scales)

library(readxl)
library(stringr)

# load data

df_to_plot <-
  read_excel("Fig6E_Data.xlsx",
             sheet = 1)

ggplot(data = df_to_plot,
       aes(x = factor(Sample,
                      levels = c("logFC_VPS45_KD",
                                 "logFC_lncRNA_KD",
                                 "logFC_C1orf54_KD",
                                 "logFC_lncRNA_KD:logFC_C1orf54_KD",
                                 "logFC_VPS45_KD:logFC_lncRNA_KD",
                                 "logFC_VPS45_KD:logFC_C1orf54_KD",
                                 "logFC_KDx3")),
           y = Coefficient,
           ymin = Coefficient - Std.Err / 2,
           ymax = Coefficient + Std.Err / 2,
           fill = 0 - log10(P))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(shape = 21,
                  size = 0.7) +
  xlab("") +
  ylim(-0.4, 0.4) +
  scale_fill_gradient(name = "-log10P",
                      low = "white",
                      # limits = c(0, 9),
                      high = "darkred") +
  # scale_colour_gr

  theme_light() +
  coord_flip()


df_to_plot <-
  df_to_plot[c(1:3, 7), ]
ggplot(data = df_to_plot,
       aes(x = factor(Sample,
                      levels = c("logFC_VPS45_KD",
                                 "logFC_lncRNA_KD",
                                 "logFC_C1orf54_KD",
                                 "logFC_KDx3")),
           y = Coefficient,
           ymin = Coefficient - Std.Err / 2,
           ymax = Coefficient + Std.Err / 2,
           fill = 0 - log10(P))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(shape = 21,
                  size = 0.7) +
  xlab("") +
  ylim(-0.4, 0.4) +
  scale_fill_gradient(name = "-log10P",
                      low = "white",
                      # limits = c(0, 9),
                      high = "darkred") +
  # scale_colour_gr

  theme_light() +
  coord_flip()
