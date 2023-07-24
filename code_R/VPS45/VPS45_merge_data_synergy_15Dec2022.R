# Siwei 15 Dec 2022
# Merge VPS45 sequencing results from single KD and 3xKD
# use sva comBat to remove batch effect
# run Kristen's synergistic model

# init
# init
library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)
library(MASS)

library(magrittr)
library(dplyr)
library(ggpubr)

library(ggplot2)
library(gplots)
library(RColorBrewer)

library(sva)
library(FactoMineR)
library(future)

library(corrplot)
library(ggrepel)
library(reshape2)
library(pheatmap)

# source functions
source("../R_synergy_analysis/code_files_Echevarria-Vargas/Melanoma_Synergy_Analysis.R")

# load data

df_raw_VPS45_single_KD <- 
  read_delim("ReadsPerGene_STAR_full.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)
# remove the A11 samples

df_VPS45_single_KD <-
  df_raw_VPS45_single_KD[, c(1, 15:ncol(df_raw_VPS45_single_KD))]

df_raw_VPS45_3xKD <-
  read_delim("VPS45_3xKD_ReadsPerGene_STAR_08Dec2022_stranded_correct.txt", 
             "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE)
# remove the A11 samples

df_VPS45_3x_KD <-
  df_raw_VPS45_3xKD[, c(1, 4:ncol(df_raw_VPS45_3xKD))]

# merge the two dataframes into one
df_master <-
  merge(x = df_VPS45_single_KD,
        y = df_VPS45_3x_KD,
        by = "Geneid")

# add gene symbols
gencode_v35_ENSG_Genename_final <- 
  read_delim("gencode.v35.ENSG.Genename.final.list", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[, c(1, 5)]
colnames(gencode_v35_ENSG_Genename_final) <-
  c("Geneid", "Gene_symbol")

gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]

df_master <-
  merge(x = df_master,
        y = gencode_v35_ENSG_Genename_final,
        by = "Geneid")
df_master <-
  df_master[!duplicated(df_master$Gene_symbol), ]

geneid_list <- df_master$Geneid
genesymbol_list <- df_master$Gene_symbol

df_master$Geneid <- NULL
df_master$Gene_symbol <- NULL

rownames(df_master) <- genesymbol_list

df_metadata <-
  data.frame(genotype = c(rep_len("C1orf54_KD", length.out = 3),
                          rep_len("EGFP", length.out = 3),
                          rep_len("lncRNA_KD", length.out = 3),
                          rep_len("VPS45_KD", length.out = 3),
                          rep_len("GG", length.out = 3),
                          rep_len("GG", length.out = 3),
                          rep_len("EGFP", length.out = 3),
                          rep_len("KDx3", length.out = 3)),
             group = factor(c(rep_len("C1orf54_KD", length.out = 3),
                              rep_len("EGFP_H12_KDx1", length.out = 3),
                              rep_len("lncRNA_KD", length.out = 3),
                              rep_len("VPS45_KD", length.out = 3),
                              rep_len("GG_H6", length.out = 3),
                              rep_len("GG_B11", length.out = 3),
                              rep_len("EGFP_H12_KDx3", length.out = 3),
                              rep_len("KDx3", length.out = 3)),
                            levels = c("EGFP_H12_KDx1",
                                       "C1orf54_KD",
                                       "lncRNA_KD", 
                                       "VPS45_KD",
                                       "KDx3",
                                       "GG_B11", 
                                       "GG_H6",
                                       "EGFP_H12_KDx3")),
             batch = c(rep_len("KDx1", length.out = 15),
                       rep_len("KDx3", length.out = 9)))

mod = model.matrix(~ as.factor(group),
                   data = df_metadata)
# registered()
# bpparam()
# 
# register(MulticoreParam(workers = 4))
# df_comBat_4_DGEList <-
#   ComBat(dat = df_master,
#          batch = as.factor(df_metadata$batch),
#          # mod = mod,
#          par.prior = T,
#          prior.plots = T,
#          BPPARAM = MulticoreParam(workers = 4,
#                                   progressbar = T))

df_comBat_4_DGEList <-
  ComBat_seq(counts = as.matrix(df_master),
             batch = as.factor(df_metadata$batch))

# remove genes that cannot be adjusted by ComBat 
# (all zero reads in at least one batch or both batches)
df_comBat_4_DGEList <- as.data.frame(df_comBat_4_DGEList)
df_comBat_4_DGEList <-
  df_comBat_4_DGEList[!(rowSums(df_comBat_4_DGEList[, 1:15] == 0) |
                         rowSums(df_comBat_4_DGEList[, 16:ncol(df_comBat_4_DGEList)] == 0)), ]

df_comBat_4_DGEList <-
  df_comBat_4_DGEList[!(rowSums(df_comBat_4_DGEList) == 0), ]

DGE_raw <-
  DGEList(counts = as.matrix(df_comBat_4_DGEList),
          samples = colnames(df_comBat_4_DGEList),
          group = df_metadata$group,
          genes = rownames(df_comBat_4_DGEList),
          remove.zeros = T)



# calc PCA use raw and corrected counts
# PCA may not look good here since we have not appplied
# the DE gene list
subsetted_df_master <-
  df_master[rownames(df_master) %in% rownames(df_comBat_4_DGEList), ]

pca_original <-
  PCA(t(subsetted_df_master),
      scale.unit = T,
      ncp = 5,
      graph = F)
fviz_pca_ind(pca_original,
             repel = T,
             habillage = as.factor(df_metadata$group),
             palette = brewer.pal(name = "Dark2", n = 8),
             show.legend = T,
             invisible = "quali",
             font.x = 10) +
  theme_light()

pca_postComBat <-
  PCA(t(df_comBat_4_DGEList),
      scale.unit = T,
      ncp = 5,
      graph = F)
fviz_pca_ind(pca_postComBat,
             repel = T,
             habillage = as.factor(df_metadata$group),
             palette = brewer.pal(name = "Dark2", n = 8),
             show.legend = T,
             invisible = "quali",
             font.x = 10) +
  theme_light()

## proceed with DEG analysis
## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_raw)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_raw <- 
  DGE_raw[((rowSums(cpm_gene_count[, df_metadata$group %in% "EGFP_H12_KDx1"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "C1orf54_KD"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "lncRNA_KD"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "VPS45_KD"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "KDx3"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "GG_B11"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "GG_H6"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "EGFP_H12_KDx3"] >= cpm_cutoff) >= 2)), ]


DGE_raw <- calcNormFactors(DGE_raw)

# use voom transformation
voom_design_matrix <- 
  model.matrix(~ 0 + as.factor(group),
               data = df_metadata)

df_voom <-
  voom(counts = DGE_raw,
       design = voom_design_matrix,
       plot = T,
       save.plot = F)

# fit, use lmFit
fit <- lmFit(df_voom,
             design = voom_design_matrix)

colnames(voom_design_matrix) <-
  str_split(string = colnames(voom_design_matrix),
            pattern = "\\)",
            simplify = T)[, 2]

cont_matrix <-
  makeContrasts(additive = (C1orf54_KD - EGFP_H12_KDx1) +
                  (lncRNA_KD - EGFP_H12_KDx1) +
                  (VPS45_KD - EGFP_H12_KDx1),
                combinational = (KDx3 - EGFP_H12_KDx1),
                synergistic = (KDx3 - C1orf54_KD - lncRNA_KD - VPS45_KD + EGFP_H12_KDx1),
                levels = voom_design_matrix)

pheatmap(t(cont_matrix),
         display_numbers = T, number_format = "%.0f",
         breaks = seq(-3, 1, by = 0.5),
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(12),
         cluster_cols = F, cluster_rows = F)

# calc coeff and STDerr
fit_contrast <-
  contrasts.fit(fit = fit,
                contrasts = cont_matrix)

# Run Empirical Bayes moderation
fit_contrast <- eBayes(fit = fit_contrast,
                       robust = T)
plotSA(fit_contrast)
summary(decideTests(fit_contrast, adjust.method = "fdr"))

# Calc DEG result table
res_list <- list()
i <- 1L
for (i in 1:length(colnames(fit_contrast$contrasts))) {
  x <- topTable(fit_contrast, 
                coef = i, 
                sort.by = "p", 
                number = Inf, 
                confint = T)
  res_list[[i]] <- x
  names(res_list)[i] <- colnames(fit_contrast$contrasts)[i]
  # write.csv(x, paste0("results/", "melanoma_analysis", "_DEGs_",
  #                     colnames(fit.cont$contrasts)[i], ".csv"))
}


# Calculate power
SE <-
  sqrt(fit_contrast$s2.post) * fit_contrast$stdev.unscaled
sig1 <- median(SE[,"additive"])
sig2 <- median(SE[,"combinational"])

# remember the function below needs to be sourced from the "Melanoma_Synergy_Analysis.R"
SE_graph <- 
  power.compare.logFC(sig1, 
                      sig2, 
                      N = 4, 
                      N_other = c(4, 6, 8, 10, 14),
                      alpha = 0.05, 
                      n_tests = 20000)
print(SE_graph)

# Determine the extent of synergy
synergy_pvalues <- res_list$synergistic$P.Value
pi1 <- 1 - qvalue(synergy_pvalues)$pi0
print(pi1)
plot.new()
text(x = 0.4, y = 0.75,
     labels = paste0("\n",round(pi1 * 100, 2),
                            "% non-null \np-values and \n",
                            round(sum(res_list$synergistic$adj.P.Val < 0.1) *
                                    100/length(res_list$synergistic$genes), 2),
                            " % of genes with \nsynergy FDR < 0.1"))
hist(synergy_pvalues)

# determine cutoff range
meanSE = mean(SE)

# adjusted names in list and fetch ensembl ID (if any), gene name, logFC, adj P value
colnames(res_list$additive)

a <- unique(res_list$additive[,c(1, 2, 8)])
b <- unique(res_list$combinational[,c(1, 2, 8)])
c <- unique(res_list$synergistic[,c(1, 2, 8)])
log2FC.matrix <- merge(x = a, 
                       y = b[!duplicated(b[, c("genes")]),],
                       all.x = TRUE, 
                       by = c("genes"))
log2FC.matrix <- merge(x = log2FC.matrix, 
                       y = c[!duplicated(c[, c("genes")]),],
                       all.x = TRUE, 
                       by = c("genes"))

colnames(log2FC.matrix) <-
  c("Gene_name",
    "Additive.logFC", "Additive.FDR",
    "Combinatorial.logFC", "Combinatorial.FDR",
    "Synergistic.logFC", "Synergistic.FDR")
log2FC.matrix <- categorize.synergy(log2FC.matrix, meanSE)
genes.per.category <- count(log2FC.matrix, "magnitude.syn")
print(genes.per.category)

# plot 
genes.per.category$category <- 
  factor(genes.per.category$magnitude.syn,
         levels = c("same", "less.up", "less.down", "more.up", "more.down"))
genes.per.category$percent <- 
  paste0(round(genes.per.category$freq *
                 100/sum(genes.per.category$freq), 0), " %")

zissou <- 
  wes_palette("Zissou1", 
              6, 
              type = "continuous")

ggplot(genes.per.category, 
       aes(x = "", 
           y = freq, 
           fill = category)) +
  geom_col() +
  coord_polar("y", 
              start = 0) +
  scale_fill_manual(values = zissou) +
  theme_void()
