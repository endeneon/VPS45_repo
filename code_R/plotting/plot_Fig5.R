# Siwei 10 Jan 2021
# plot Fig 5, S8-

# init
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(readxl)

# plot 5b_mfr

df_plot <- read_excel("plot_data.xlsx",
                      sheet = "Fig5b_mfr")

ggplot(df_plot,
       aes(x = Day,
           y = MFR,
           fill = Genotype)) +
  geom_boxplot(outlier.alpha = 0,
               position = "dodge2") +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              # position = "identity",
  #              position = position_dodge(0.75),
  #              dotsize = 0.6) +
  # geom_point()
  # geom_jitter(aes(colour = Genotype),
  #             stat = "identity",
  #             position = "jitter",
  #             # width = 0.1,
  #             size = 1) +
  # scale_colour_manual(values = c("darkred", "black")) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  # ylim(0, 12) +
  xlab("Days post NGN2 induction") +
  ylab("Mean firing rate (Hz)") +
  theme_classic() #+
  # theme(legend.position = "none")

# t test for AA vs GG

t.test(x = df_plot$MFR[(df_plot$Genotype == "AA") &
                         (df_plot$Day == "Day 50")],
       y = df_plot$MFR[(df_plot$Genotype == "GG") &
                         (df_plot$Day == "Day 50")],
       alternative = "t",
       paired = F,
       var.equal = F)

t.test(x = df_plot$MFR[(df_plot$Genotype == "AA") &
                         (df_plot$Day == "Day 53")],
       y = df_plot$MFR[(df_plot$Genotype == "GG") &
                         (df_plot$Day == "Day 53")],
       alternative = "t",
       paired = F,
       var.equal = F)

t.test(x = df_plot$MFR[(df_plot$Genotype == "AA") &
                         (df_plot$Day == "Day 56")],
       y = df_plot$MFR[(df_plot$Genotype == "GG") &
                         (df_plot$Day == "Day 56")],
       alternative = "t",
       paired = F,
       var.equal = F)


# plot 5b_NOB

df_plot <- read_excel("plot_data.xlsx",
                      sheet = "Fig5b_nob")

ggplot(df_plot,
       aes(x = Day,
           y = NOB,
           fill = Genotype)) +
  geom_boxplot(outlier.alpha = 0,
               position = "dodge2") +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               # position = "identity",
               position = position_dodge(0.75),
               dotsize = 0.6) +
  # geom_point()
  # geom_jitter(aes(colour = Genotype),
  #             stat = "identity",
  #             position = "jitter",
  #             # width = 0.1,
  #             size = 1) +
  # scale_colour_manual(values = c("darkred", "black")) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  ylim(0, 3000) +
  xlab("Days post NGN2 induction") +
  ylab("Number of bursts") +
  theme_classic() +
  theme(legend.position = "none")

# t test for AA vs GG

t.test(x = df_plot$NOB[(df_plot$Genotype == "AA") &
                         (df_plot$Day == "Day 50")],
       y = df_plot$NOB[(df_plot$Genotype == "GG") &
                         (df_plot$Day == "Day 50")],
       alternative = "t",
       paired = F,
       var.equal = F)

t.test(x = df_plot$NOB[(df_plot$Genotype == "AA") &
                         (df_plot$Day == "Day 53")],
       y = df_plot$NOB[(df_plot$Genotype == "GG") &
                         (df_plot$Day == "Day 53")],
       alternative = "t",
       paired = F,
       var.equal = F)

t.test(x = df_plot$NOB[(df_plot$Genotype == "AA") &
                         (df_plot$Day == "Day 56")],
       y = df_plot$NOB[(df_plot$Genotype == "GG") &
                         (df_plot$Day == "Day 56")],
       alternative = "t",
       paired = F,
       var.equal = F)
