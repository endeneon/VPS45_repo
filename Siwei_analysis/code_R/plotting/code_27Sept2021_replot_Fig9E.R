# Siwei 13 Sept 2021
# Siwei 27 Sept 2021
# replot Fig S9E using another layout
# make bubble plot for Fig S9E


# init
library(readr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

library(readxl)

library(stringr)

# load data
## C1orf54
df_raw_data <- read_excel("FigS10D_data.xlsx",
                          sheet = 1)
df_raw <- df_raw_data
# df_raw <- df_raw_data[1:40, ]
# df_raw <- df_raw_data[41:80, ]

df_raw$`-log10FDR` <- 0 - log10(df_raw$FDR)

# KD_condition_raw_order <-
#   unique(as.factordf_raw$KD_Condition))
df_raw$KD_Condition <- as.factor(df_raw$KD_Condition)
# df_raw$Rank <- rep_len(c(1:10),
#                        length.out = 40)
# df_raw$Rank <- as.factor(df_raw$Rank)
df_raw$Part <- as.factor(df_raw$Part)

# df_raw <- df_raw[str_detect(df_raw$description,
#                             pattern = "neur"), ]

df_raw <- df_raw[order(df_raw$enrichmentRatio,
                       decreasing = T), ]
# df_raw <- df_raw[order(df_raw$`-log10FDR`,
#                        decreasing = T), ]
df_raw <- df_raw[1:40, ]

ggplot(df_raw, aes(x = description,
                   y = KD_Condition,
                   size = enrichmentRatio,
                   fill = `-log10FDR`)) +
  scale_fill_gradient(low = "white",
                      # limits = c(0, 7),
                      high = "darkred") +
  scale_radius(limits = c(0, max(df_raw_data$enrichmentRatio))) +
  geom_point(shape = 21) +
  # geom_text(aes(label = description),
  #           hjust = 0,
  #           vjust = 0,
  #           size = 3) +
  xlab("C1orf54 Knockdown") +
  # ylab("Rank") +
  # scale_y_discrete(limits = rev(levels(df_raw$Rank))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = unit(c(1, 3, 1, 1),
                           units = "cm"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.title = element_text(size = 10)) #+

## VPS45KD
df_raw_data <- read_excel("FigS10D_data.xlsx",
                          sheet = 2)
df_raw <- df_raw_data
# df_raw <- df_raw_data[1:40, ]
# df_raw <- df_raw_data[41:80, ]

df_raw$`-log10FDR` <- 0 - log10(df_raw$FDR)

# KD_condition_raw_order <-
#   unique(as.factordf_raw$KD_Condition))
df_raw$KD_Condition <- as.factor(df_raw$KD_Condition)
# df_raw$Rank <- rep_len(c(1:10),
#                        length.out = 40)
# df_raw$Rank <- as.factor(df_raw$Rank)
df_raw$Part <- as.factor(df_raw$Part)

# df_raw <- df_raw[str_detect(df_raw$description,
#                             pattern = "neur"), ]

df_raw <- df_raw[order(df_raw$enrichmentRatio,
                       decreasing = T), ]
# df_raw <- df_raw[order(df_raw$`-log10FDR`,
#                        decreasing = T), ]
df_raw <- df_raw[1:40, ]

ggplot(df_raw, aes(x = description,
                   y = KD_Condition,
                   size = enrichmentRatio,
                   fill = `-log10FDR`)) +
  scale_fill_gradient(low = "white",
                      # limits = c(0, 7),
                      high = "darkred") +
  scale_radius(limits = c(0, max(df_raw_data$enrichmentRatio))) +
  geom_point(shape = 21) +
  # geom_text(aes(label = description),
  #           hjust = 0,
  #           vjust = 0,
  #           size = 3) +
  xlab("VPS45 Knockdown") +
  # ylab("Rank") +
  # scale_y_discrete(limits = rev(levels(df_raw$Rank))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = unit(c(1, 3, 1, 1),
                           units = "cm"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.title = element_text(size = 10)) #+

## lncRNAKD
df_raw_data <- read_excel("FigS10D_data.xlsx",
                          sheet = 3)
df_raw <- df_raw_data
# df_raw <- df_raw_data[1:40, ]
# df_raw <- df_raw_data[41:80, ]

df_raw$`-log10FDR` <- 0 - log10(df_raw$FDR)

# KD_condition_raw_order <-
#   unique(as.factordf_raw$KD_Condition))
df_raw$KD_Condition <- as.factor(df_raw$KD_Condition)
# df_raw$Rank <- rep_len(c(1:10),
#                        length.out = 40)
# df_raw$Rank <- as.factor(df_raw$Rank)
df_raw$Part <- as.factor(df_raw$Part)

# df_raw <- df_raw[str_detect(df_raw$description,
#                             pattern = "neur"), ]

df_raw <- df_raw[order(df_raw$enrichmentRatio,
                       decreasing = T), ]
# df_raw <- df_raw[order(df_raw$`-log10FDR`,
#                        decreasing = T), ]
df_raw <- df_raw[1:40, ]

ggplot(df_raw, aes(x = description,
                   y = KD_Condition,
                   size = enrichmentRatio,
                   fill = `-log10FDR`)) +
  scale_fill_gradient(low = "white",
                      # limits = c(0, 7),
                      high = "darkred") +
  scale_radius(limits = c(0, max(df_raw_data$enrichmentRatio))) +
  geom_point(shape = 21) +
  # geom_text(aes(label = description),
  #           hjust = 0,
  #           vjust = 0,
  #           size = 3) +
  xlab("lncRNA Knockdown") +
  # ylab("Rank") +
  # scale_y_discrete(limits = rev(levels(df_raw$Rank))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = unit(c(1, 3, 1, 1),
                           units = "cm"),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.title = element_text(size = 10)) #+
