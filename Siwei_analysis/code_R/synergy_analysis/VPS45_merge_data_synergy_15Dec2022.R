# Siwei 15 Dec 2022
# Merge VPS45 sequencing results from single KD and 3xKD
# use sva comBat to remove batch effect
# run Kristen's synergistic model

# init
# init
# library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)
library(MASS)

# library(magrittr)
# library(dplyr)
# library(ggpubr)
#
# library(ggplot2)
# library(gplots)
# library(RColorBrewer)
#
library(sva)
library(FactoMineR)
library(future)
#
# library(corrplot)
# library(ggrepel)
# library(reshape2)
# library(pheatmap)

pacman::p_load(limma, edgeR, pheatmap, RColorBrewer, colorspace,
               ggplot2, ggpubr, qvalue,  plyr, wesanderson,
               GSEABase, grid, scales, WebGestaltR, stringr)

# source functions
source(file = "Synergy_functions_only.R")

# load data

df_raw_VPS45_single_KD <-
  read_delim("VPS45/ReadsPerGene_STAR_full.txt",
             "\t", escape_double = FALSE, trim_ws = TRUE)
# remove the A11 samples

df_VPS45_single_KD <-
  df_raw_VPS45_single_KD[, c(1, 15:ncol(df_raw_VPS45_single_KD))]

df_raw_VPS45_3xKD <-
  read_delim("VPS45/VPS45_3xKD_ReadsPerGene_STAR_08Dec2022_stranded_correct.txt",
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

df_master$Geneid <-
  str_split(string = df_master$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_master <-
  df_master[!duplicated(df_master$Geneid), ]

rm(list = c("df_raw_VPS45_3xKD", "df_raw_VPS45_single_KD",
            "df_VPS45_3x_KD", "df_VPS45_single_KD"))

geneid_list <- df_master$Geneid
df_master$Geneid <- NULL
rownames(df_master) <- geneid_list

#
# rownames(df_master) <- gene_name_list
# # add gene symbols
# gencode_v35_ENSG_Genename_final <-
#   read_delim("gencode.v35.ENSG.Genename.final.list",
#              delim = "\t", escape_double = FALSE,
#              col_names = FALSE, trim_ws = TRUE)
# gencode_v35_ENSG_Genename_final <-
#   gencode_v35_ENSG_Genename_final[, c(1, 5)]
# colnames(gencode_v35_ENSG_Genename_final) <-
#   c("Geneid", "Gene_symbol")
#
# gencode_v35_ENSG_Genename_final <-
#   gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
# gencode_v35_ENSG_Genename_final <-
#   gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]
#
# df_master <-
#   merge(x = df_master,
#         y = gencode_v35_ENSG_Genename_final,
#         by = "Geneid")
# df_master <-
#   df_master[!duplicated(df_master$Gene_symbol), ]
#
# geneid_list <- df_master$Geneid
# gene_name_list <- df_master$Gene_symbol
#
# df_master$Geneid <- NULL
# df_master$Gene_symbol <- NULL
#
# rownames(df_master) <- gene_name_list

df_metadata_raw <-
  data.frame(samples = colnames(df_master),
             genotype = c(rep_len("C1orf54_KD", length.out = 3),
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

df_metadata <-
  df_metadata_raw
#
df_comBat_4_DGEList <- df_master
#
# df_metadata <-
#   df_metadata[match(colnames(df_comBat_4_DGEList),
#                     df_metadata$samples), ]
## !!! the combinatorial effects of GG or 3xKD have to be separately evaluated
## choose the corresponding columns from df_master and df_metadata_raw here
## Remember to set different analysis names when choosing different sample combinations
##



#### Select 1xKD vs 3xKD group, use EGFP from 1xKD as control,
#### KDx3 as combinatorial, use EGFP_H12_KDx3 as the control for KDx3
experiment.title <- "KDx1_vs_KDx3"

df_comBat_4_DGEList <-
  df_master[, !((df_metadata_raw$group %in% "GG_B11") |
                  (df_metadata_raw$group %in% "GG_H6"))]

df_metadata <-
  df_metadata_raw[!((df_metadata_raw$group %in% "GG_B11") |
                      (df_metadata_raw$group %in% "GG_H6")), ]
df_metadata <-
  df_metadata[match(colnames(df_comBat_4_DGEList),
                    df_metadata$samples), ]

## !! H12_sg_2 since it not group well !!
#### remove H12_sg_2 for better contrast
df_comBat_4_DGEList <-
  df_comBat_4_DGEList[, !(df_metadata$samples %in% "H12_sg_2_")]
df_metadata <-
  df_metadata[!(df_metadata$samples %in% "H12_sg_2_"), ]
df_metadata <-
  df_metadata[match(colnames(df_comBat_4_DGEList),
                    df_metadata$samples), ]

# df_comBat_4_DGEList <-
#   df_comBat_4_DGEList[, c(-8, -14)]
# df_metadata <-
#   df_metadata[c(-8, -14), ]
# rownames(df_metadata) <- colnames(df_comBat_4_DGEList)

##!! Need to redefine the level of df_metadata$group to exclude levels with zero element
df_metadata$group <-
  factor(df_metadata$group,
         levels = c("EGFP_H12_KDx1",
                    "C1orf54_KD",
                    "lncRNA_KD",
                    "VPS45_KD",
                    "EGFP_H12_KDx3",
                    "KDx3"))

df_metadata$group <-
  factor(df_metadata$group,
         levels = c("EGFP_H12_KDx1",
                    "C1orf54_KD",
                    "lncRNA_KD",
                    "VPS45_KD",
                    "GG_B11",
                    "GG_H6",
                    "EGFP_H12_KDx3",
                    "KDx3"))


######## FINISHED 1xKD vs 3xKD value assignment

#### Select 1xKD vs GG group, use EGFP from 1xKD as control, GG as combinatorial
# experiment.title <- "KDx1_vs_GG"
#### PENDING

### run ComBat_seq
### Note: do not use the ComBat() since it introduces negative values, it should
### be used on normalised reads (e.g., FPKM, TPM)
# df_comBat_4_DGEList_backup <- df_comBat_4_DGEList
# df_comBat_4_DGEList <- df_comBat_4_DGEList_backup
#
df_comBat_4_DGEList <-
  ComBat_seq(counts = as.matrix(df_comBat_4_DGEList),
             batch = as.factor(df_metadata$batch))

# remove genes that cannot be adjusted by ComBat
# (all zero reads in at least one batch or both batches)
df_comBat_4_DGEList <- as.data.frame(df_comBat_4_DGEList)
df_comBat_4_DGEList <-
  df_comBat_4_DGEList[!(rowSums(df_comBat_4_DGEList[, 1:12]) == 0 |
                         rowSums(df_comBat_4_DGEList[, 13:ncol(df_comBat_4_DGEList)]) == 0), ]

df_comBat_4_DGEList <-
  df_comBat_4_DGEList[!(rowSums(df_comBat_4_DGEList[, 1:18]) == 0 |
                          rowSums(df_comBat_4_DGEList[, 19:ncol(df_comBat_4_DGEList)]) == 0), ]


df_comBat_4_DGEList <-
  df_comBat_4_DGEList[!(rowSums(df_comBat_4_DGEList) == 0), ]


DGE_raw <-
  DGEList(counts = as.matrix(df_comBat_4_DGEList),
          samples = colnames(df_comBat_4_DGEList),
          group = df_metadata$group,
          genes = rownames(df_comBat_4_DGEList),
          remove.zeros = T)

anno <- read.csv("data_code_brennand_nat_protocol/anno.csv")
rownames(anno) <- anno$ensembl
anno <-
  anno[match(rownames(DGE_raw),
             rownames(anno)), ]
DGE_raw$genes <- anno

## proceed with DEG analysis
## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_raw)

cpm_gene_count <-
  cpm_gene_count[rownames(cpm_gene_count) %in% geneid_1267_list, ]

## If analysing synergy-driving gene exp, may skip this part,
## continue at calculating CPM

# calc PCA use raw and corrected counts
# PCA may not look good here since we have not appplied
# the DE gene list
# subsetted_df_master <-
#   df_master[rownames(df_master) %in% rownames(df_comBat_4_DGEList), ]
#
# pca_original <-
#   PCA(t(subsetted_df_master),
#       scale.unit = T,
#       ncp = 5,
#       graph = F)
# fviz_pca_ind(pca_original,
#              repel = T,
#              habillage = as.factor(df_metadata$group),
#              palette = brewer.pal(name = "Dark2", n = 8),
#              show.legend = T,
#              invisible = "quali",
#              font.x = 10) +
#   theme_light()
#
# pca_postComBat <-
#   PCA(t(cpm_gene_count),
#       scale.unit = T,
#       ncp = 5,
#       graph = F)

pca_postComBat <-
  prcomp(x = t(cpm_gene_count),
         center = T,
         scale. = T)

fviz_pca_ind(pca_postComBat,
             repel = T,
             habillage = as.factor(df_metadata$group),
             palette = brewer.pal(name = "Dark2", n = 8),
             show.legend = T,
             invisible = "quali",
             font.x = 10) +
  theme_light()



## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 0.4

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_raw <-
  DGE_raw[((rowSums(cpm_gene_count[, df_metadata$group %in% "EGFP_H12_KDx1"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "C1orf54_KD"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "lncRNA_KD"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "VPS45_KD"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "KDx3"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "GG_B11"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "GG_H6"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, df_metadata$group %in% "EGFP_H12_KDx3"] >= cpm_cutoff) >= 3)), ]


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

## !! the elements in count matrix needs to change as the group name changes
## (e.g., KDx3 changed to GG)
cont_matrix <-
  makeContrasts(additive = (C1orf54_KD + lncRNA_KD + VPS45_KD - 3 * EGFP_H12_KDx1),
                combinatorial = (KDx3 - EGFP_H12_KDx3),
                synergistic = (KDx3 -
                                 C1orf54_KD - lncRNA_KD - VPS45_KD -
                                 EGFP_H12_KDx3 +
                                 3 * EGFP_H12_KDx1),
                levels = voom_design_matrix)

# pheatmap(t(cont_matrix),
#          display_numbers = T, number_format = "%.0f",
#          breaks = seq(-3, 1, by = 0.5),
#          color = colorRampPalette(rev(brewer.pal(n = 10, name = "diverge_hcl")))(12),
#          cluster_cols = F, cluster_rows = F)

pheatmap(t(cont_matrix),
         display_numbers = T, number_format = "%.0f",
         breaks = seq(-3, 1, by = 0.5),
         color = diverge_hcl(n = 12,
                             h = c(255, 330),
                             # c = 70,
                             l = c(40, 90)),
         number_color = "black",
         cluster_cols = F, cluster_rows = F)

# calc coeff and STDerr
fit_contrast <-
  contrasts.fit(fit = fit,
                contrasts = cont_matrix)

# Run Empirical Bayes moderation
fit_contrast <- eBayes(fit = fit_contrast)

plotSA(fit_contrast)
summary(decideTests(fit_contrast, adjust.method = "fdr"))
# additive combinatorial synergistic
# Down        869          2025          57
# NotSig    20308         17493       20961
# Up          220          1879         379


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
  write.csv(x, paste0("results/", "VPS45_analysis", "_DEGs_",
                      colnames(fit_contrast$contrasts)[i], ".csv"))
}


# Calculate power
SE <-
  sqrt(fit_contrast$s2.post) * fit_contrast$stdev.unscaled
colnames(SE)
sig1 <- median(SE[,"additive"])
sig2 <- median(SE[,"combinatorial"])

# remember the function below needs to be sourced from the "Melanoma_Synergy_Analysis.R"
SE_graph <-
  power.compare.logFC(sig1,
                      sig2,
                      N = 4,
                      N_other = c(4, 6, 8, 10, 12, 14),
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
                                    100/length(res_list$synergistic$Gene_name), 2),
                            " % of genes with \nsynergy FDR < 0.1"))
hist(synergy_pvalues)

# determine cutoff range
meanSE = mean(SE)

# adjusted names in list and fetch ensembl ID (if any), gene name, logFC, adj P value
colnames(res_list$additive)

a <- unique(res_list$additive[,c(1, 2, 4, 9)])
b <- unique(res_list$combinatorial[,c(2, 4, 9)])
c <- unique(res_list$synergistic[,c(2, 4, 9)])
log2FC.matrix <- merge(x = a,
                       y = b[!duplicated(b[, c("Gene_name")]),],
                       all.x = TRUE,
                       by = c("Gene_name"))
log2FC.matrix <- merge(x = log2FC.matrix,
                       y = c[!duplicated(c[, c("Gene_name")]),],
                       all.x = TRUE,
                       by = c("Gene_name"))

colnames(log2FC.matrix) <-
  c("Ensembl", "Gene_name",
    "Additive.logFC", "Additive.FDR",
    "Combinatorial.logFC", "Combinatorial.FDR",
    "Synergistic.logFC", "Synergistic.FDR")

# remove incorrectly formatted gene names
log2FC.matrix <- log2FC.matrix[-c(1:25), ]
log2FC.matrix <- log2FC.matrix[!is.na(log2FC.matrix$Gene_name), ]

log2FC.matrix <-
  log2FC.matrix[order(log2FC.matrix$Gene_name), ]
write.table(x = log2FC.matrix,
            file = "log2FC_matrix.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

log2FC.matrix <- categorize.synergy(log2FC.matrix, meanSE)
genes.per.category <-
  plyr::count(log2FC.matrix,
              vars = "magnitude.syn")
print(genes.per.category)
# magnitude.syn  freq  category percent
# 1     less.down  2597 less.down    14 %
# 2       less.up  2585   less.up    14 %
# 3     more.down  1288 more.down     7 %
# 4       more.up  1171   more.up     6 %
# 5          same 10858      same    59 %

# plot
genes.per.category$category <-
  factor(genes.per.category$magnitude.syn,
         levels = c("same", "less.up", "less.down", "more.up", "more.down"))
genes.per.category$percent <-
  paste0(round(genes.per.category$freq *
                 100/sum(genes.per.category$freq), 0), " %")

# zissou <-
#   wes_palette("Zissou1",
#               6,
#               type = "continuous")
zissou <-
  brewer.pal(name = "Set1",
              n = 6)

# plot the pie chart of genes per category
ggplot(genes.per.category,
       aes(x = "",
           y = freq,
           fill = category)) +
  geom_col() +
  coord_polar("y",
              start = 0) +
  scale_fill_manual(values = c("#7FC97F",
                               # "#F69679","#F69679",
                               "#F69679", "#8781BD",
                               "#EE0000", "#0088FF")) +
  theme_void() #+
  margin(1, 1, 1, 1, unit = "in")

# plot heatmaps of the log2FC in the additive and
# the combinatorial comparisons for each synergy category.
for (i in 1:length(levels(log2FC.matrix$magnitude.syn))) {
  breaks <- c(seq(-6, -0.3, by = 0.1),
              seq(0.3, 6, by = 0.1))
  breaks <- append(breaks, -9, 0)
  breaks <- append(breaks, 9)
  tmp <- log2FC.matrix[log2FC.matrix$magnitude.syn ==
                         levels(log2FC.matrix$magnitude.syn)[i],
                       c("Additive.logFC","Combinatorial.logFC")]
  h <-
  pheatmap(tmp,
           kmeans_k = 30,
           cellwidth = 70, cellheight = 5,
           border_color = NA,
           breaks=breaks,
           cluster_cols = F,
           show_rownames = F,
           color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(117),
           main = paste0("logFC expected vs. measured:\n",
                         levels(log2FC.matrix$magnitude.syn)[i]))
  print(h)
}


#### Enrichment analysis #### Schrode & Brennand et al Nature Protocols

## All comparisons: GSEA ##


## 23.	Create a list containing the gene set groups of interest
# in the required format using ids2indices(), geneIds() and getGmt().

# In this example these are manually curated gene set groups,
# saved in the “genesets” folder:
# disorder.gmt, behavior.gmt, connectivity.gmt, head.gmt, neural.gmt,
# postsynapse.gmt, presynapse.gmt and synapse.gmt.

gs.list <-
  list("disorder" = ids2indices(geneIds(getGmt("genesets/disorder.gmt")),
                                identifiers = df_voom$genes$Gene_name),
       "behavior" = ids2indices(geneIds(getGmt("genesets/behavior.gmt")),
                                identifiers = df_voom$genes$Gene_name),
       "connectivity" = ids2indices(geneIds(getGmt("genesets/connectivity.gmt")),
                                    identifiers = df_voom$genes$Gene_name),
       "head" = ids2indices(geneIds(getGmt("genesets/head.gmt")),
                            identifiers = df_voom$genes$Gene_name),
       "neural" = ids2indices(geneIds(getGmt("genesets/neural.gmt")),
                              identifiers = df_voom$genes$Gene_name),
       "postsynapse" = ids2indices(geneIds(getGmt("genesets/postsynapse.gmt")),
                                   identifiers = df_voom$genes$Gene_name),
       "presynapse" = ids2indices(geneIds(getGmt("genesets/presynapse.gmt")),
                                  identifiers = df_voom$genes$Gene_name),
       "synapse" = ids2indices(geneIds(getGmt("genesets/synapse.gmt")),
                               identifiers = df_voom$genes$Gene_name))

# Create a custom color palette with at least as many colors as gene set groups to be tested.
# Use the Dark2 palette
catcols <-
  brewer.pal(n = 8, name = "Accent")
names(catcols) <-
  c("behaviour", "disorder", "connectivity", "head",
    "neural", "postsynapse", "presynapse", "synapse")

show_col(catcols)

##	Run gene set enrichment for all comparisons using camera.

# Loop through all contrasts in the cont.matrix object.
# Perform enrichment and visualize as scatter
# and bar plots using the custom cameraplusplots() function.

j <- 1
camera.res.list <- list()
for (j in 1:length(colnames(cont_matrix))) {
  print(paste0("Contrast: ", colnames(cont_matrix)[j]))
  # pdf(paste0("results/", experiment.title, "-10_GSEA-", colnames(cont_matrix)[j], "-plots.pdf"))
  camera.res <-
    cameraplusplots(contrast = cont_matrix[ , j],
                    genesetlist = gs.list,
                    vobject = df_voom,
                    design = voom_design_matrix,
                    catcolors = catcols,
                    title = paste0(colnames(cont_matrix)[j]))
  # dev.off()
  camera.res.list[[j]] <- camera.res
  names(camera.res.list)[j] <- colnames(cont_matrix)[j]
  write.csv(data.frame(camera.res), paste0("results/", experiment.title,
                                           "_GSEA-", colnames(cont_matrix)[j], ".csv"))
}

# Plot legend.
# pdf(paste0("results/", experiment.title, "-10_GSEA-plot-legend.pdf"))
plot(1,type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
legend("center", names(catcols), cex = 1.2, fill = catcols)
# dev.off()


## Specific gene subsets: ORA ##

## 26. Create synergistic gene subsets to be analyzed.

# Adjust FDR cutoff depending on the subtlety of the synergistic effect.
# In this example a cutoff of synergistic FDR < 1% was chosen.

log2FC.matrix.sub <- subset(log2FC.matrix, Synergistic.FDR < 0.05)
colnames(log2FC.matrix.sub) <-
  c("Gene_name", "Ensembl",
    "Additive.logFC", "Additive.FDR",
    "Combinatorial.logFC", "Combinatorial.FDR",
    "Synergistic.logFC", "Synergistic.FDR",
    "magnitude.syn")

#------------
# Further possibly interesting subsets (commented out):

### Synergy genes with combinatorial FDR < 5%
#log2FC.matrix.sub <- subset(log2FC.matrix, Combinatorial.FDR < 0.05)
### Genes with synergistic Fold Change > 2 or < 0.5
#log2FC.matrix.sub <- subset(log2FC.matrix, Synergistic.logFC > 1.5 | Synergistic.logFC < -1.5)
#------------


## 27. Stratify chosen subset by synergy category using the custom stratify.by.syn.cat() function.

syn.cat.list <- stratify.by.syn.cat(log2FC.matrix.sub)


## 28.	Define reference genes as all genes analyzed.
# Since lowly expressed genes were filtered out at the beginning of the protocol,
# this can be interpreted as all expressed genes.

allgenes <- as.character(DGE_raw$genes$Gene_name)


## 29.	Create a list of paths for the gene sets of interest.
# WebgestaltR, which is used for over-representation analysis,
# requires gene sets to be provided in a different format than camera.

gs.list <- list(
  "disorder" = "genesets/disorder.gmt",
  "behavior" = "genesets/behavior.gmt",
  "connectivity" = "genesets/connectivity.gmt",
  "head" = "genesets/head.gmt",
  "neural" = "genesets/neural.gmt",
  "postsynapse" = "genesets/postsynapse.gmt",
  "presynapse" = "genesets/presynapse.gmt",
  "synapse" = "genesets/synapse.gmt")


## 30.	Run over-representation analysis for “more up” and “more down” synergy categories using the WebGestaltR() function.

# Loop through the “more.up” and “more.down” vectors in the previously created list object syn.cat.list.
# Create a data frame for all results.
# Visualize using the custom oraplot() function.

## more.down
for (i in 3:3) {
  tryCatch({
    ora <- data.frame(matrix(ncol = 11, nrow = 0))
    ora.list <- list()
    goi <- syn.cat.list[[i]]
    for (j in 1:length(gs.list)) {
      tryCatch({
        ora.s <- WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             interestGene = goi,
                             interestGeneType = "genesymbol",
                             referenceGene = allgenes,
                             referenceGeneType = "genesymbol",
                             enrichDatabase = "others",
                             enrichDatabaseFile = file.path(gs.list[j]),
                             enrichDatabaseType = "genesymbol",
                             sigMethod = "top", topThr = 50, minNum = 3,
                             isOutput = F)
        ora.list[[j]] <- ora.s
        names(ora.list)[j] <- names(gs.list)[j]
        ora.list[[j]]$category <- names(ora.list[j])
        colnames(ora) <- names(ora.list[[1]])
        ora <- rbind.data.frame(ora, ora.list[[j]])
      },
      error = function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })
    }
    # write.csv(ora, paste0("results/", experiment.title, "_ORA-",
    #                       levels(log2FC.matrix$magnitude.syn)[i], ".csv"))

    g <- oraplot(ora, catcols,
                 paste0(levels(log2FC.matrix$magnitude.syn)[i]))
    # pdf(paste0("results/", experiment.title, "-11_ORA-",
    #            levels(log2FC.matrix$magnitude.syn)[i], "-plots.pdf"))
    print(g)
    # dev.off()
  },
  error = function(e){
    cat("ERROR :",conditionMessage(e), "\n")
  })
}

## make gitter plot of the more up/down regulated genes
# syn.cat.list[[3]]: more.down
# syn.cat.list[[4]]: more.up
ora$dirNeglogFDR <-
  0 - log10(ora$FDR)

ggplot(aes(x = category,
           y = dirNeglogFDR,
           color = category),
       data = ora) +
  scale_color_manual(values = catcols) +
  geom_jitter(aes(size = size,
                  alpha = dirNeglogFDR),
              pch = 19,
              show.legend = F) +
  scale_size_continuous(range = c(1, 8)) +
  scale_alpha_continuous(range = c(0.4, 1)) +
  geom_hline(yintercept = 1.3,
             color = "red",
             alpha = 1) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(0, 15),
                     oob = squish,
                     labels = abs) +
  labs(x = "Gene set categories", y = "-log10(FDR)",
       title = "more.down") +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


## more.up
for (i in 4:4) {
  tryCatch({
    ora <- data.frame(matrix(ncol = 11, nrow = 0))
    ora.list <- list()
    goi <- syn.cat.list[[i]]
    for (j in 1:length(gs.list)) {
      tryCatch({
        ora.s <- WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             interestGene = goi,
                             interestGeneType = "genesymbol",
                             referenceGene = allgenes,
                             referenceGeneType = "genesymbol",
                             enrichDatabase = "others",
                             enrichDatabaseFile = file.path(gs.list[j]),
                             enrichDatabaseType = "genesymbol",
                             sigMethod = "top", topThr = 50, minNum = 3,
                             isOutput = F)
        ora.list[[j]] <- ora.s
        names(ora.list)[j] <- names(gs.list)[j]
        ora.list[[j]]$category <- names(ora.list[j])
        colnames(ora) <- names(ora.list[[1]])
        ora <- rbind.data.frame(ora, ora.list[[j]])
      },
      error = function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })
    }
    # write.csv(ora, paste0("results/", experiment.title, "_ORA-",
    #                       levels(log2FC.matrix$magnitude.syn)[i], ".csv"))

    g <- oraplot(ora, catcols,
                 paste0(levels(log2FC.matrix$magnitude.syn)[i]))
    # pdf(paste0("results/", experiment.title, "-11_ORA-",
    #            levels(log2FC.matrix$magnitude.syn)[i], "-plots.pdf"))
    print(g)
    # dev.off()
  },
  error = function(e){
    cat("ERROR :",conditionMessage(e), "\n")
  })
}

## make gitter plot of the more up/down regulated genes
# syn.cat.list[[3]]: more.down
# syn.cat.list[[4]]: more.up
ora$dirNeglogFDR <-
  0 - log10(ora$FDR)

ggplot(aes(x = category,
           y = dirNeglogFDR,
           color = category),
       data = ora) +
  scale_color_manual(values = catcols) +
  geom_jitter(aes(size = size,
                  alpha = dirNeglogFDR),
              pch = 19,
              show.legend = F) +
  scale_size_continuous(range = c(1, 8)) +
  scale_alpha_continuous(range = c(0.4, 1)) +
  geom_hline(yintercept = 1.3,
             color = "red",
             alpha = 1) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(0, 15),
                     oob = squish,
                     labels = abs) +
  labs(x = "Gene set categories", y = "-log10(FDR)",
       title = "more.up") +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

save.image("VPS45_synergy_data_set_22Dec2022.RData")
