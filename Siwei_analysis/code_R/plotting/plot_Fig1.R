# Siwei 10 Jan 2021
# plot Fig 1, S1-S2

# init
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(readxl)

# plot Fig1b

df_plot <- read_excel("plot_data.xlsx",
                      sheet = "Fig1b")

ggplot(df_plot,
       aes(x = Ref_ratio...10,
           y = Ref_ratio...17)) +
  geom_point(size = 1) +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  xlim(0, 1) +
  ylim(0, 1) +
  ylab("Reference allele ratio in iN-Glut") +
  xlab("Reference allele ration in NGN2 neuron") +
  stat_smooth(method = "lm",
              se = F,
              fullrange = T) +
  ggtitle("Pearson's R = 0.809, P = 5.75x10-5") +
  theme_classic()

cor.test(df_plot$Ref_ratio...10,
         df_plot$Ref_ratio...17,
         method = "p",
         alternative = "t")

# plot Fig S2a
df_plot <- read_excel("plot_data.xlsx",
                      sheet = "Fig_S2a")
df_plot$Ref_ratio <- df_plot$REF_C / (df_plot$REF_C + df_plot$ALT_C)
df_plot$neglogP <- 0 - log10(df_plot$pVal)
df_plot$FDR_threshold <- "FDR >= 0.05"
df_plot$FDR_threshold[df_plot$FDR < 0.05] <- "FDR < 0.05"
df_plot$FDR_threshold <- factor(df_plot$FDR_threshold)

ggplot(df_plot,
       aes(x = Ref_ratio,
           y = neglogP,
           colour = FDR_threshold)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = c("red", "black")) +
  xlim(0, 1) +
  ylim(0, 50) +
  xlab("Ratio of the reference allele") +
  ylab("-log10 P value") +
  theme_classic() +
  theme(legend.position = "none")

# plot Fig S2b
df_plot <- read_excel("plot_data.xlsx",
                      sheet = "Fig_S2b")

ggplot(df_plot,
       aes(x = Ref_ratio...10,
           y = Ref_ratio...17)) +
  geom_point(size = 1) +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  xlim(0, 1) +
  ylim(0, 1) +
  ylab("Reference allele ratio in iN-Glut") +
  xlab("Reference allele ration in NGN2 neuron") +
  stat_smooth(method = "lm",
              se = F,
              fullrange = T) +
  ggtitle("Pearson's R = 0.410, P = 4.23x10-14") +
  theme_classic()

cor.test(df_plot$Ref_ratio...10,
         df_plot$Ref_ratio...17,
         method = "p",
         alternative = "t")
