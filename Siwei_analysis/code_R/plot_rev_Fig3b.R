# Siwei 23 Dec 2020
# make revised Fig 3b
# use both AA/GG and AG/GG

# init
library(ggplot2)
library(scales)

#
Fig3b_data_ref <- Fig3b_data_ref[, 1]
Fig3b_data <- merge(Fig3b_data_ref, Fig3b_data,
                    by = "Geneid")


Fig3b_data$colour <- "darkred"
Fig3b_data$colour <- "dodgerblue3"

Fig3b_data$alpha <-
  0 - log10(Fig3b_data$FDR)
Fig3b_data$alpha <-
  Fig3b_data$alpha / max(Fig3b_data$alpha)
Fig3b_data$SE <-
  abs(Fig3b_data$logFC / (-0.862 + (sqrt(0.743 - 2.404*log(Fig3b_data$PValue)))))




ggplot(Fig3b_data, aes(x = START,
                       y = 0 - logFC)) +
  geom_point(colour = Fig3b_data$colour,
             fill = Fig3b_data$colour,
             size = 1,
             shape = 23) +
  geom_errorbar(aes(x = START,
                    ymin = 0 - (logFC - SE / 2),
                    ymax = 0 - (logFC + SE / 2)),
                colour = Fig3b_data$colour,
                width = 0.5,
                position = "dodge2",
                size = 0.5) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  scale_x_continuous(labels = scales::scientific_format(digits = 4)) +
  xlab("Chr1 position") +
  ylab("Expression fold change") +
  ylim(-1.2, 1.2) +
  scale_alpha_continuous(name = "-log10FDR") +
  theme_classic()
