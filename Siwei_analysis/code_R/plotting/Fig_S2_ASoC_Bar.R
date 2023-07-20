# Siwei 31 Jan 2022
# re-plot Fig S2 ASoC bar plot

# init
library(ggplot2)
library(RColorBrewer)

library(readxl)
library(scales)

S2_ASoC_bar <- read_excel("S2_ASoC_bar.xlsx")

ggplot(S2_ASoC_bar,
       aes(x = eQTL,
           y = Percentage,
           fill = Category)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))
