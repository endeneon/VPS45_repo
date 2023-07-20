# Siwei 01 Sept 2021
# Re-plot Fig 6G

# init
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(colorspace)

library(scales)

library(readxl)
library(stringr)

# load data
## neuro
df_to_plot <-
  read_excel("Fig6G_Data.xlsx",
             sheet = 1)

ggplot(df_to_plot,
       aes(x = Diff,
           y = Value,
           fill = Trend)) +
  geom_bar(position = "fill",
           stat = "identity") +
  geom_text(aes(label = Value),
            position = position_stack(vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("") +
  ylab("Percentage") +

  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  facet_grid(. ~ factor(Category,
                        levels = c("lncRNA-KD",
                                   "VPS45-KD",
                                   "C1ORF54-KD")),
             switch = "x")

## synaptic
df_to_plot <-
  read_excel("Fig6G_Data.xlsx",
             sheet = 2)

ggplot(df_to_plot,
       aes(x = Diff,
           y = Value,
           fill = Trend)) +
  geom_bar(position = "fill",
           stat = "identity") +
  geom_text(aes(label = Value),
            position = position_stack(vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("") +
  ylab("Percentage") +

  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0)) +
  facet_grid(. ~ factor(Category,
                        levels = c("lncRNA-KD",
                                   "VPS45-KD",
                                   "C1ORF54-KD")),
             switch = "x")
dev.off()

