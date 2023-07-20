# Siwei 27 Aug 2021
# Re-plot Fig S8 with merged data

# init
library(ggplot2)
library(RColorBrewer)

library(readxl)

# > RColorBrewer::brewer.pal(3, "Dark2")
# [1] "#1B9E77" "#D95F02" "#7570B3"

# load data
# weighted mean firing rate
df_to_plot <-
  read_excel("Fig_S8_rev_data.xlsx",
             sheet = "Weighted_MFR")

ggplot(data = df_to_plot,
       aes(x = Day,
           colour = Genotype,
           y = Value)) +
  geom_line(size = 1) +
  geom_point(shape = 21,
             fill = "white") +
  scale_color_discrete(type = brewer.pal(3, "Dark2")) +
  ylab("Hz") +
  theme_classic() +
  ggtitle("Weighted Mean Firing Rate (Hz)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

# No_of_bursts
df_to_plot <-
  read_excel("Fig_S8_rev_data.xlsx",
             sheet = 2)

ggplot(data = df_to_plot,
       aes(x = Day,
           colour = Genotype,
           y = Value)) +
  geom_line(size = 1) +
  geom_point(shape = 21,
             fill = "white") +
  scale_color_discrete(type = brewer.pal(3, "Dark2")) +
  ylab("# of events/10 min") +
  theme_classic() +
  ggtitle("Number of Bursts") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

# Burst Duration
df_to_plot <-
  read_excel("Fig_S8_rev_data.xlsx",
             sheet = 3)

ggplot(data = df_to_plot,
       aes(x = Day,
           colour = Genotype,
           y = Value)) +
  geom_line(size = 1) +
  geom_point(shape = 21,
             fill = "white") +
  scale_color_discrete(type = brewer.pal(3, "Dark2")) +
  ylab("Duration (s)") +
  theme_classic() +
  ggtitle("Burst Duration") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

# Network Burst Frequency
df_to_plot <-
  read_excel("Fig_S8_rev_data.xlsx",
             sheet = 4)

ggplot(data = df_to_plot,
       aes(x = Day,
           colour = Genotype,
           y = Value)) +
  geom_line(size = 1) +
  geom_point(shape = 21,
             fill = "white") +
  scale_color_discrete(type = brewer.pal(3, "Dark2")) +
  ylab("Hz") +
  theme_classic() +
  ggtitle("Network Burst Frequency") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

# Network Burst Duration
df_to_plot <-
  read_excel("Fig_S8_rev_data.xlsx",
             sheet = 5)

ggplot(data = df_to_plot,
       aes(x = Day,
           colour = Genotype,
           y = Value)) +
  geom_line(size = 1) +
  geom_point(shape = 21,
             fill = "white") +
  scale_color_discrete(type = brewer.pal(3, "Dark2")) +
  ylab("Duration (s)") +
  theme_classic() +
  ggtitle("Network Burst Duration") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

# Synchrony Index
df_to_plot <-
  read_excel("Fig_S8_rev_data.xlsx",
             sheet = 6)

ggplot(data = df_to_plot,
       aes(x = Day,
           colour = Genotype,
           y = Value)) +
  geom_line(size = 1) +
  geom_point(shape = 21,
             fill = "white") +
  scale_color_discrete(type = brewer.pal(3, "Dark2")) +
  ylab("Duration (s)") +
  theme_classic() +
  ggtitle("Index") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))
