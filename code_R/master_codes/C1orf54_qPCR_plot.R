# Siwei 29 Mar 2021
# make C1orf54 qPCR plot

# init
library(readr)
library(ggplot2)

# 
qPCR_table <- read_delim("qPCR_table.txt", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

df_to_plot <- as.data.frame(t(qPCR_table),
                            stringsAsFactors = F)
colnames(df_to_plot) <- df_to_plot[1, ]
df_to_plot <- df_to_plot[-1, ]

df_to_plot$C1orf54_qPCR_frag <- as.numeric(df_to_plot$C1orf54_qPCR_frag)
df_to_plot$Assigned <- as.numeric(df_to_plot$Assigned)

df_to_plot$pseudo_cpm <- 
  df_to_plot$C1orf54_qPCR_frag * 1e6 / df_to_plot$Assigned

df_to_plot$cell_line <- c(rep_len("A11", length.out = 6),
                          rep_len("H12", length.out = 6))
df_to_plot$treatment <- c(rep_len("C1orf54KD", length.out = 3),
                          rep_len("EGFP", length.out = 3),
                          rep_len("C1orf54KD", length.out = 3),
                          rep_len("EGFP", length.out = 3))

df_to_plot$cell_line <- factor(df_to_plot$cell_line)
df_to_plot$treatment <- factor(df_to_plot$treatment)


ggplot(df_to_plot, 
       aes(x = treatment,
           y = pseudo_cpm,
           # colour = cell_line,
           fill = cell_line)) +
  geom_dotplot(binaxis = "y",
               stackdir = "center") +
  geom_boxplot(position = "dodge2") +
  ylab("pseudo_cpm_qPCR_frag") +
  theme_classic()


ggplot(df_to_plot, 
       aes(# x = treatment,
           x = pseudo_cpm,
           colour = cell_line,
           fill = cell_line)) +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center") +
  geom_bar(position = "dodge") +

  theme_classic()
