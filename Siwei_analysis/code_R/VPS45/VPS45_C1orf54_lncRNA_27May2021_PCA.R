# Siwei 27 May 2021
# edgeR analysis for all VPS45-related RNASeq data
# make PCA and hierachical clustering plots


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

##### Not using edgeR built-in GO analysis #####
## generate gene expression list only

gene_count_table <- 
  read_delim("ReadsPerGene_STAR_RNASeq_C1orf54_VPS45_lncRNA_GG.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
# gene_count_table$X17 <- NULL
# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

gene_list <- gene_count_table$Geneid
gene_count_table$Geneid <- NULL
gene_count_table <- gene_count_table[, -c(25:27)] ## remove GG columns 
rownames(gene_count_table) <- gene_list

## assemble gencode table
gencode_v35_ENSG_Genename_final <- 
  read_delim("gencode.v35.ENSG.Genename.final.list", 
             "\t", escape_double = FALSE, col_names = FALSE, 
             trim_ws = TRUE)

colnames(gencode_v35_ENSG_Genename_final) <- 
  c("Geneid", "CHR", "START", "END", "Gene_symbol")
gencode_v35_ENSG_Genename_final$Geneid <- 
  str_split(gencode_v35_ENSG_Genename_final$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]



## make the raw DGEList from the count matrix
DGE_Raw <- DGEList(as.matrix(gene_count_table),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45",
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

DGE_cpm <- cpm(DGE_Raw)
DGE_cpm <- as.data.frame(DGE_cpm, 
                         stringsAsFactors = F)
# DGE_cpm[rownames(DGE_cpm) %in% "11311", ]

samples <- data.frame(genotype = c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45",
                                   "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                      cell.line = c("A11", "A11", "A11", "A11", "A11", "A11",
                                    "A11", "A11", "A11", "A11", "A11", "A11",
                                    "H12", "H12", "H12", "H12", "H12", "H12",
                                    "H12", "H12", "H12", "H12", "H12", "H12"))

# design_samples <- model.matrix(~ factor(samples$cell.line) +
#                                  factor(samples$genotype))

design_samples <- model.matrix(~ factor(samples$genotype,
                                        levels = c("AA_EGFP", "AA_C1orf54",
                                                   "AA_lncRNA", "AA_VPS45")))

cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_lncRNA"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_C1orf54"] >= cpm_cutoff) >= 6)
          , ]

DGE_Raw <- calcNormFactors(DGE_Raw,
                           method = "TMM")#c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw <- estimateDisp(DGE_Raw, 
                        design = design_samples,
                        robust = T)

## make MDS and BCV plots to observe the clustering pattern of samples
## with the cell-type element regressed out as blocking factor.
plotMDS(DGE_Raw)
plotBCV(DGE_Raw)

MDS_data <- plotMDS(DGE_Raw)




DGE_QLF <- glmQLFit(DGE_Raw,
                    design = design_samples,
                    robust = T)
DGE_QLF <- glmQLFTest(DGE_QLF,
                      coef = 3) # refer to the structure of the design matrix
DGE_QLF <- glmQLFTest(DGE_QLF,
                      contrast = c(-1, 0.33, 0.33, 0.33))
DGE_results <- DGE_QLF$table
DGE_results$FDR <- p.adjust(DGE_results$PValue,
                            method = "fdr")
hist(DGE_results$FDR, 
     breaks = 1000, xlim = c(0, 0.1))
DGE_results <- 
  DGE_results[order(DGE_results$PValue), ]
top_gene_list <-
  rownames(DGE_results)[1:(as.integer(nrow(DGE_results) / 5))]

save(top_gene_list,
     file = "VPS45_top_gene_list_4_PCA_plot.RData")
# DGE_QLF <- glmQLFTest(DGE_QLF,
#                       contrast = c(-1, 0, 1)) # refer to the structure of the design matrix
table(decideTests(DGE_QLF))

### load full gene table
gene_count_table_full <- 
  read_delim("ReadsPerGene_STAR_RNASeq_C1orf54_VPS45_lncRNA_GG.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)
# gene_count_table$X17 <- NULL
# rename and rearrange Geneid
gene_count_table_full$Geneid <- str_split(gene_count_table_full$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table_full <- 
  gene_count_table_full[!duplicated(gene_count_table_full$Geneid), ]

gene_list <- gene_count_table_full$Geneid
gene_count_table_full$Geneid <- NULL
rownames(gene_count_table_full) <- gene_list

### subset genes that are only included in the top_gene_list
gene_count_table_full <-
  gene_count_table_full[rownames(gene_count_table_full) %in% top_gene_list, ]

DGE_full <- DGEList(as.matrix(gene_count_table_full),
                   genes = rownames(gene_count_table_full),
                   samples = colnames(gene_count_table_full),
                   remove.zeros = T)

### calculate cpm for PCA
cpm_4_pca <- as.data.frame(cpm(DGE_full))
cpm_4_pca <- as.data.frame(t(cpm_4_pca))

pca_cpm <- prcomp(x = cpm_4_pca,
                  center = T,
                  scale. = T,
                  tol = 0)

fviz_pca_ind(pca_cpm,
             col.ind = "contrib",
             gradient.cols = c("darkblue", "grey50", "red"),
             repel = T)

pca_plot_data <- fviz_pca_ind(pca_cpm,
                              col.ind = "contrib",
                              gradient.cols = c("darkblue", "grey50", "red"),
                              repel = T)

df_pca_plot <- data.frame(x = pca_plot_data$data$x,
                          y = pca_plot_data$data$y,
                          name = pca_plot_data$data$name,
                          stringsAsFactors = F)
df_pca_plot$genotype = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                "AA_VPS45", "AA_VPS45", "AA_VPS45",
                                "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                "AA_VPS45", "AA_VPS45", "AA_VPS45",
                                "GG_EGFP", "GG_EGFP", "GG_EGFP"),
                              levels = c("AA_EGFP", "AA_C1orf54",
                                         "AA_lncRNA", "AA_VPS45",
                                         "GG_EGFP"))

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
                                "#7570B3")) +
  xlab("PC1 = 32.4%") +
  ylab("PC2 = 22.1%") +
  theme_classic() +
  geom_text_repel()

### make hierachical clustering

heatmap.2(as.matrix(gene_count_table),
          dendrogram = "column",
          col = "cm.colors",
          scale = "row",
          trace = "none",
          labRow = "")
  
