# plot Fig 2J

# init
library(readxl)
library(ggplot2)

library(RColorBrewer)

library(stringr)

# load data

# down_regulated genes
df_plot <-
  read_excel("Fig_2J_table.xlsx",
             sheet = "AA_down")

df_plot <-
  read_excel("Fig2J_new_17Apr2023.xlsx")

df_plot <-
  df_plot[1:10, ]

df_plot$order <-
  rev(1:nrow(df_plot))
df_plot$GO_Terms <-
  str_c(df_plot$geneSet,
        df_plot$description,
        sep = "-")
df_plot$GO_Terms <-
  factor(df_plot$GO_Terms,
         levels = df_plot$GO_Terms)



ggplot(df_plot,
       aes(y = GO_Terms,
           x = 0 - log10(FDR),
           fill = enrichmentRatio)) +
  geom_point(shape = 23,
             size = 2) +
  geom_vline(xintercept = (0 - log10(0.05)),
             colour = "red",
             alpha = 0.5) +
  xlim(0, max(0 - log10(df_plot$FDR))) +
  xlab("-log10(FDR)") +
  scale_y_discrete(limits = rev(df_plot$GO_Terms)) +
  scale_fill_gradient2(low = "white",
                       high = "darkred",
                       limits = c(0, 35),
                       name = "Enrichment Ratio") +
  theme_classic() #+
  ggtitle("AA_down")


# up_regulated genes
df_plot <-
  read_excel("Fig_2J_table.xlsx",
             sheet = "AA_up")

df_plot$order <-
  rev(1:nrow(df_plot))
df_plot$GO_Terms <-
  str_c(df_plot$geneSet,
        df_plot$description,
        sep = "-")
df_plot$GO_Terms <-
  factor(df_plot$GO_Terms,
         levels = df_plot$GO_Terms)

ggplot(df_plot,
       aes(y = GO_Terms,
           x = 0 - log10(FDR),
           fill = enrichmentRatio)) +
  geom_point(shape = 23,
             size = 2) +
  geom_vline(xintercept = (0 - log10(0.05)),
             colour = "red",
             alpha = 0.5) +
  xlim(0, max(0 - log10(df_plot$FDR))) +
  scale_y_discrete(limits = rev(df_plot$GO_Terms)) +
  scale_fill_gradient2(low = "white",
                       high = "darkred",
                       limits = c(0, max(df_plot$enrichmentRatio)),
                       name = "") +
  theme_classic() +
  ggtitle("AA_up")
