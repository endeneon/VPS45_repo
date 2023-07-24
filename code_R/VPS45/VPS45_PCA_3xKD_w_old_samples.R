# Siwei 08 Dec 2022

# make PCA plot using anchor gene list from 27 May 2021



# init
library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(ggrepel)

# load top gene list
load("VPS45_top_gene_list_4_PCA_plot.RData")


### load May 2021 gene table
gene_count_table_May2021 <- 
  read_delim("ReadsPerGene_STAR_RNASeq_C1orf54_VPS45_lncRNA_GG.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)
# gene_count_table$X17 <- NULL
# rename and rearrange Geneid
gene_count_table_May2021$Geneid <- str_split(gene_count_table_May2021$Geneid,
                                          pattern = "\\.",
                                          simplify = T)[, 1]
gene_count_table_May2021 <- 
  gene_count_table_May2021[!duplicated(gene_count_table_May2021$Geneid), ]

gene_list <- gene_count_table_May2021$Geneid
# gene_count_table_May2021$Geneid <- NULL
rownames(gene_count_table_May2021) <- gene_list

# load Dec 2022 gene table
gene_count_table_Dec2022 <- 
  read_delim("VPS45_3xKD_ReadsPerGene_STAR_08Dec2022_stranded_correct.txt", 
             "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE)
gene_count_table_Dec2022$`A11-EGFP-1_` <- NULL
gene_count_table_Dec2022$`A11-EGFP-2_` <- NULL
gene_count_table_Dec2022$Geneid <- 
  str_split(gene_count_table_Dec2022$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
gene_count_table_Dec2022 <- 
  gene_count_table_Dec2022[!duplicated(gene_count_table_Dec2022$Geneid), ]

gene_list <- gene_count_table_Dec2022$Geneid
# gene_count_table_Dec2022$Geneid <- NULL
rownames(gene_count_table_Dec2022) <- gene_list


gene_count_table_full <-
  merge(x = gene_count_table_May2021,
        y = gene_count_table_Dec2022,
        by = "Geneid")
gene_list <- gene_count_table_full$Geneid
gene_count_table_full$Geneid <- NULL
rownames(gene_count_table_full) <- gene_list

### subset genes that are only included in the top_gene_list
gene_count_table_full <-
  gene_count_table_full[rownames(gene_count_table_full) %in% top_gene_list, ]

colnames(gene_count_table_full) <-
  str_replace_all(string = colnames(gene_count_table_full),
                  pattern = "_$",
                  replacement = "")

DGE_full <- DGEList(as.matrix(gene_count_table_full),
                    genes = rownames(gene_count_table_full),
                    samples = colnames(gene_count_table_full),
                    remove.zeros = T)

### calculate cpm for PCA
cpm_4_pca <- as.data.frame(cpm(DGE_full))
cpm_4_pca <- as.data.frame(t(cpm_4_pca))

pca_cpm <- 
  prcomp(x = cpm_4_pca,
         center = T,
         scale. = T,
         tol = 0)

pca_plot_data <- 
  fviz_pca_ind(pca_cpm,
               col.ind = "contrib",
               gradient.cols = c("darkblue", "grey50", "red"),
               repel = T)

df_pca_plot <- 
  data.frame(x = pca_plot_data$data$x,
             y = pca_plot_data$data$y,
             name = pca_plot_data$data$name,
             stringsAsFactors = F)


df_pca_plot$genotype <- 
  factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
           "AA_EGFP", "AA_EGFP", "AA_EGFP",
           "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
           "AA_VPS45", "AA_VPS45", "AA_VPS45",
           "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
           "AA_EGFP", "AA_EGFP", "AA_EGFP",
           "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
           "AA_VPS45", "AA_VPS45", "AA_VPS45",
           "GG_EGFP", "GG_EGFP", "GG_EGFP",
           "GG_EGFP_Dec2022", "GG_EGFP_Dec2022", "GG_EGFP_Dec2022", 
           "AA_EGFP_Dec2022", "AA_EGFP_Dec2022", "AA_EGFP_Dec2022", 
           "AA_3xKD_Dec2022", "AA_3xKD_Dec2022", "AA_3xKD_Dec2022"),
         levels = c("AA_EGFP", "AA_C1orf54",
                    "AA_lncRNA", "AA_VPS45",
                    "GG_EGFP",
                    "AA_EGFP_Dec2022", "AA_3xKD_Dec2022", "GG_EGFP_Dec2022"))


### make PCA plot
ggplot(df_pca_plot, aes(x = x,
                        y = y,
                        label = name,
                        colour = genotype)) +
  geom_point() +
  scale_color_manual(name = "Genotypes and\n Knockdown Genes",
                     values = c("#1B9E77",
                                "#E7298A",
                                "#D95F02",
                                "#666666",
                                "#7570B3",
                                "#660000",
                                "#006600",
                                "#000066")) +
  xlab("PC1 = 46.7%") +
  ylab("PC2 = 16.1%") +
  theme_classic() +
  geom_text_repel()
