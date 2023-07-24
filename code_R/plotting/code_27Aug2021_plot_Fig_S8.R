# Siwei 27 Aug 2021
# Re-plot Fig S8 with averaged out numbers

# init
library(ggplot2)
library(RColorBrewer)

library(readxl)

# > RColorBrewer::brewer.pal(3, "Dark2")
# [1] "#1B9E77" "#D95F02" "#7570B3"

# load data
# weighted mean firing rate
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = "Weighted_MFR")

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("Hz") +
  theme_classic() +
  ggtitle("Weighted Mean Firing Rate (Hz)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))


# number of bursts
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = 2)

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("# of events in 10 min") +
  theme_classic() +
  ggtitle("Number of Bursts") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))



# burst duration
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = 3)

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("Duration (s)") +
  theme_classic() +
  ggtitle("Burst Duration") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))


# Number of spikes per burst
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = 4)

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("# of spikes") +
  theme_classic() +
  ggtitle("Number of Spikes per Burst") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))


# network burst duration
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = 5)

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("Duration (s)") +
  xlim(35, 74) +
  theme_classic() +
  ggtitle("Network Burst Duration") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

# Number of spikes per network burst
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = 6)

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("# of spikes") +
  xlim(35, 74) +
  theme_classic() +
  ggtitle("Number of Spikes per Network Burst") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))


# synchrony index
df_to_plot <-
  read_excel("Fig_S8_data.xlsx",
             sheet = 7)

# remove all even rows
df_to_plot <-
  df_to_plot[seq(1, nrow(df_to_plot), 2), ]

ggplot(df_to_plot,
       aes(x = Days)) +
  geom_line(aes(y = AA),
            color = "#1B9E77",
            size = 1) +
  geom_line(aes(y = AG),
            color = "#D95F02",
            size = 1) +
  geom_line(aes(y = GG),
            color = "#7570B3",
            size = 1) +
  geom_point(aes(x = Days,
                 y = AA),
             color = "#1B9E77",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = AG),
             color = "#D95F02",
             size = 1,
             shape = 21,
             fill = "white") +
  geom_point(aes(x = Days,
                 y = GG),
             color = "#7570B3",
             size = 1,
             shape = 21,
             fill = "white") +
  ylab("Index") +
  # xlim(35, 74) +
  theme_classic() +
  ggtitle("Synchrony Index") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))
