# Siwei 14 Feb 2023
# make log10CPM expression bar plot of genes

# init
library(readxl)
library(ggplot2)
library(RColorBrewer)

## NGN2_neuron_ATAC_seq
df_raw <-
  read_excel("Neuron_iPSC-Specific marker RNA-seq_expression_NGN2 neurons.xlsx",
             sheet = 3)
gene_list <- df_raw$Gene

df_raw$Gene <- NULL
rownames(df_raw) <-
  gene_list

df_to_plot <-
  data.frame(log10CPM = unlist(df_raw),
             Gene = rep(x = unlist(gene_list),
                        times = ncol(df_raw)))

ggplot(df_to_plot,
       aes(x = factor(Gene,
                      levels = gene_list),
           y = log10CPM,
           fill = factor(Gene))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(width = 0.1,
              size = 1) +
  ylab("log2CPM") +
  # ylim(0, max(df_to_plot$log10CPM) * 1.2) +
  # xlab("") +
  scale_fill_manual(values = brewer.pal(n = 9,
                                        name = "Set1"),
                    name = "Gene") +
  ylim(-5, 12) +
  theme_classic() +
  ggtitle("NGN2_neuron_from_ATAC_seq") +
  theme(axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 315,
                                   vjust = 1,
                                   hjust = 0))


## iPS
df_raw <-
  read_excel("Neuron_iPSC-Specific marker RNA-seq_expression_NGN2 neurons.xlsx",
             sheet = 4)
gene_list <- df_raw$Gene

df_raw$Gene <- NULL
rownames(df_raw) <-
  gene_list

df_to_plot <-
  data.frame(log10CPM = unlist(df_raw),
             Gene = rep(x = unlist(gene_list),
                        times = ncol(df_raw)))

ggplot(df_to_plot,
       aes(x = factor(Gene,
                      levels = gene_list),
           y = log10CPM,
           fill = as.factor(Gene))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(width = 0.1,
              size = 1) +
  ylab("log2CPM") +
  # ylim(0, max(df_to_plot$log10CPM) * 1.2) +
  # xlab("") +
  ylim(-5, 12) +
  scale_fill_manual(values = brewer.pal(n = 9,
                                        name = "Set1"),
                    name = "Gene") +
  # scale_colour_manual(values = brewer.pal(n = 6,
  #                                         name = "Set2"),
  #                     name = "Clone") +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("iPS Cells") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 315,
                                   vjust = 1,
                                   hjust = 0))


## NGN2-CRISPR
df_raw <-
  read_excel("Neuron_iPSC-Specific marker RNA-seq_expression_NGN2 neurons.xlsx",
             sheet = 5)
gene_list <- df_raw$Gene

df_raw$Gene <- NULL
rownames(df_raw) <-
  gene_list

df_to_plot <-
  data.frame(log10CPM = unlist(df_raw),
             Gene = rep(x = unlist(gene_list),
                        times = ncol(df_raw)))

ggplot(df_to_plot,
       aes(x = factor(Gene,
                      levels = gene_list),
           y = log10CPM,
           fill = factor(Gene))) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(width = 0.1,
              size = 1) +
  ylab("log2CPM") +
  # ylim(0, max(df_to_plot$log10CPM) * 1.2) +
  # xlab("") +
  ylim(-5, 12) +
  scale_fill_manual(values = brewer.pal(n = 9,
                                        name = "Set1"),
                    name = "Gene") +
  # scale_colour_manual(values = brewer.pal(n = 6,
  #                                         name = "Set2"),
  #                     name = "Clone") +
  # position_dodge2(width = 2) +
  theme_classic() +
  ggtitle("NGN2_neuron_SNP_editing") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 315,
                                   vjust = 1,
                                   hjust = 0))
