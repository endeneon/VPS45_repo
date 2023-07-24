# Siwei 19 Sept 2021
# Plot Fig S1 purity graph

# init
library(ggplot2)
library(RColorBrewer)
library(colorspace)

library(scales)

library(readxl)

## load data
df_to_plot <- read_excel("Fig_S1_purity.xlsx")
df_to_plot$Line <- as.factor(df_to_plot$Line)

ggplot(df_to_plot,
       aes(x = Line,
           y = Purity)) +
  geom_bar(position = "dodge",
           stat = "identity",
           fill = "darkblue") +
  xlab("Cell line") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()
