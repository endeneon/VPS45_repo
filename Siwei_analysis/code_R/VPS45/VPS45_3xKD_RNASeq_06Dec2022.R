# Siwei 06 Dec 2022
# edgeR analysis for VPS45 3xKD RNASeq data

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

library(corrplot)
library(ggrepel)
library(reshape2)

# # load data
# gene_count_table <- 
#   read_delim("VPS45_3xKD_ReadsPerGene_STAR_06Dec2022.txt", 
#              "\t", 
#              escape_double = FALSE, 
#              trim_ws = TRUE)
# 
# gene_count_table <- 
#   read_delim("VPS45_3xKD_ReadsPerGene_STAR_07Dec2022.txt", 
#              "\t", 
#              escape_double = FALSE, 
#              trim_ws = TRUE)

gene_count_table <- 
  read_delim("VPS45_3xKD_ReadsPerGene_STAR_08Dec2022_stranded_correct.txt", 
             "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE)

gencode_v35_ENSG_Genename_final <- 
  read_delim("gencode.v35.ENSG.Genename.final.list", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[, c(1, 5)]
colnames(gencode_v35_ENSG_Genename_final) <-
  c("Geneid", "Gene_symbol")

# gencode_v35_ENSG_Genename_final$Geneid <-
#   str_split(gencode_v35_ENSG_Genename_final$Geneid,
#                pattern = "\\.",
#                simplify = T)[, 1]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]

df_raw <-
  merge(x = gene_count_table,
        y = gencode_v35_ENSG_Genename_final,
        by = "Geneid")

# df_raw <- gene_count_table

# df_raw <- 
#   df_raw[!duplicated(df_raw$Geneid), ]
# rownames(df_raw) <- df_raw$Geneid

df_raw <-
  df_raw[!duplicated(df_raw$Gene_symbol), ]
rownames(df_raw) <- df_raw$Gene_symbol

# rownames(df_raw) <- df_raw$Geneid

df_4_DGE <- df_raw
df_4_DGE$Geneid <- NULL
df_4_DGE$Gene_symbol <- NULL

rownames(df_4_DGE) <- df_raw$Gene_symbol

df_4_DGE$`A11-EGFP-1_` <- NULL
df_4_DGE$`A11-EGFP-2_` <- NULL
# df_4_DGE$

df_4_DGE <-
  df_4_DGE[, c(4, 5, 6, 
               7, 8, 9, 
               1, 2, 3)]

# df_4_DGE$H12_EGF_1_ <- NULL # H12_EGF_1 showed as an outlier on MDS

## make metadata df
df_metadata <-
  data.frame(samples = colnames(df_4_DGE))
df_metadata$genotype <-
  factor(x = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
               "AA_3xKD", "AA_3xKD", "AA_3xKD",
               "GG_EGFP", "GG_EGFP", "GG_EGFP"),
         levels = c("AA_EGFP", "AA_3xKD", "GG_EGFP"))
df_metadata$cell_line <-
  factor(x = c(rep_len("H12",
                       length.out = 6),
               rep_len("B11",
                       length.out = 3)),
         levels = c("H12", "B11"))

# df_metadata$genotype <-
#   factor(x = c("AA_EGFP", "AA_EGFP", 
#                "AA_3xKD", "AA_3xKD", "AA_3xKD",
#                "GG_EGFP", "GG_EGFP", "GG_EGFP"),
#          levels = c("AA_EGFP", "AA_3xKD", "GG_EGFP"))

### make DGEList
df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = df_metadata$samples,
          group = df_metadata$genotype,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE),
                stringsAsFactors = F)
# hist(rowmeans(as.matrix(cpm_gene_count)),
#      breaks = 10000)


cpm_cutoff <- 1 # 17431 genes left

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
df_DGE <- 
  df_DGE[(rowSums(cpm_gene_count[, df_metadata$genotype %in% 
                                   "GG_EGFP"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count[, df_metadata$genotype %in% 
                                      "AA_EGFP"] >= cpm_cutoff) >= 2) |
               (rowSums(cpm_gene_count[, df_metadata$genotype %in% 
                                         "AA_3xKD"] >= cpm_cutoff) >= 2) 
          , ]

retained_gene_list <- rownames(df_DGE)

df_DGE <- calcNormFactors(df_DGE)
plotMDS(df_DGE)

cpm_gene_count_filtered <-
  cpm_gene_count[rownames(cpm_gene_count) %in% 
                   rownames(df_DGE)
                 , ]

# cpm_gene_count[rownames(cpm_gene_count) %in% "ENSG00000118292.9", ]
# cpm_gene_count[rownames(cpm_gene_count) %in% "ENSG00000285154", ]
# 
# cpm_gene_count[rownames(cpm_gene_count) %in% "ENSG00000285184.3", ] # lncRNA
# cpm_gene_count[rownames(cpm_gene_count) %in% "ENSG00000136631.15", ] # VPS45
# cpm_gene_count[rownames(cpm_gene_count) %in% "ENSG00000118292.9", ] # C1orf54
# 
# df_raw$Geneid[str_detect(string = df_raw$Geneid,
#             pattern = "ENSG00000285184")]
# df_raw$Geneid[str_detect(string = df_raw$Geneid,
#                          pattern = "ENSG00000136631")]
# df_raw$Geneid[str_detect(string = df_raw$Geneid,
#                          pattern = "ENSG00000118292")]
# 
# gencode_v35_ENSG_Genename_final[gencode_v35_ENSG_Genename_final$X1 %in% "ENSG00000118292.9", ]

df_design_matrix <-
  model.matrix(~ genotype,
               data = df_metadata)

df_DGE <-
  estimateDisp(df_DGE,
               design = df_design_matrix,
               robust = T)
plotBCV(df_DGE)

# find genes different between ANY of the three groups
# (maximise group isolation)
result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)
result_QLM <-
  glmQLFTest(result_QLM,
             coef = 2:3)
summary(decideTestsDGE(result_QLM))
topTags(result_QLM)

results_table <-
  result_QLM$table
results_table$FDR <-
  p.adjust(p = results_table$PValue,
           method = "fdr")
top_genes_list <-
  rownames(results_table)[results_table$FDR < 0.05]

# AA_3xKD vs AA_EGFP
result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)
result_QLM <-
  glmQLFTest(result_QLM,
             coef = 2)
summary(decideTestsDGE(result_QLM))
topTags(result_QLM)

# GG_EGFP vs AA_EGFP
result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)
result_QLM <-
  glmQLFTest(result_QLM,
             coef = 3)
summary(decideTestsDGE(result_QLM))
topTags(result_QLM)


#### make PCA plot
df_4_PCA <-
  as.data.frame(result_QLM$fitted.values)

df_4_PCA <-
  cpm_gene_count
df_4_PCA <-
  df_4_PCA[rownames(df_4_PCA) %in% top_genes_list, ]

df_4_PCA <-
  cpm_gene_count_filtered

df_4_PCA <-
  df_DGE$counts

df_4_PCA <-
  ComBat_seq(counts = as.matrix(df_4_PCA),
             batch = df_metadata$cell_line)
# df_4_PCA <-
#   df_4_PCA

pca_res_adj <- 
  PCA(t(df_4_PCA), 
      scale.unit = T, 
      ncp = 10, 
      graph = F)


pca_res_adj <-
  prcomp(x = t(df_4_PCA),
         retx = T,
         center = T,
         scale. = T)

fviz_pca_ind(pca_res_adj,
             repel = T,
             habillage = df_metadata$genotype,
             palette = brewer.pal(name = "Set1", n = 3),
             show.legend = F,
             invisible = "quali", font.x = 10) +
  theme_light() 

# calc MDS

df_4_PCA <-
  cpm_gene_count_filtered

df_4_PCA <-
  t(df_4_PCA)

df_MDS <-
  df_4_PCA %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()

df_MDS <-
  df_4_PCA %>%
  dist() %>%
  isoMDS(maxit = 200) %>%
  .$points %>%
  as_tibble()

df_MDS <-
  df_4_PCA %>%
  dist() %>%
  sammon(niter = 200) %>%
  .$points %>%
  as_tibble()

colnames(df_MDS) <-
  c("Dim_1", "Dim_2")

ggscatter(df_MDS, 
          x = "Dim_1", y = "Dim_2", 
          label = rownames(df_4_PCA),
          size = 1,
          repel = TRUE)



#########
# re_analyse the df matrix by filtering the df_4_DGE using
# the retained_gene_list list, attach the manually counted 
# confirmed_C1orf54_transcript

df_4_DGE <-
  df_4_DGE[rownames(df_4_DGE) %in% retained_gene_list, ]
# > colnames(df_4_DGE)
# [1] "H12_EGF_1_" "H12_EGF_2_" "H12_EGF_3_" "H12_sg_1_"  "H12_sg_2_"  "H12_sg_3_" 
# [7] "B11_EGF_1_" "B11_EGF_2_" "B11_EGF_3_"

### !!! after MATCH !!!
df_4_DGE <-
  as.data.frame(rbind(df_4_DGE,
                      c(2, 16, 3, 1, 0, 3, 1, 6, 2)))
rownames(df_4_DGE)[nrow(df_4_DGE)] <- 
  "confirmed_C1orf54_transcript"

# df_4_DGE <-
#   as.data.frame(rbind(df_4_DGE,
#                       c(232, 274, 255, 
#                         133, 145, 158,
#                         256, 220, 219)))
# rownames(df_4_DGE)[nrow(df_4_DGE)] <- 
#   "VPS45_ENST00000369130"

df_4_DGE <-
  as.data.frame(rbind(df_4_DGE,
                      c(401, 293, 243, 
                        143, 150, 92,
                        205, 163, 216)))
rownames(df_4_DGE)[nrow(df_4_DGE)] <- 
  "VPS45_ENSE00001071807"


df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = df_metadata$samples,
          group = df_metadata$genotype,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE),
                stringsAsFactors = F)

gencode_v35_ENSG_Genename_final[str_detect(string = gencode_v35_ENSG_Genename_final$Geneid,
                                           pattern = "ENSG00000285184")
                                , ] # AC244033.2


df_to_plot <-
  data.frame(samples = colnames(cpm_gene_count),
             CPM = unlist(cpm_gene_count[rownames(cpm_gene_count) %in% 
                                           "AC244033.2"
                                         , ]),
             genotype = df_metadata$genotype,
             stringsAsFactors = F)

df_to_plot <-
  data.frame(samples = colnames(cpm_gene_count),
             CPM = unlist(cpm_gene_count[rownames(cpm_gene_count) %in% 
                                           "confirmed_C1orf54_transcript"
                                         , ]),
             genotype = df_metadata$genotype,
             stringsAsFactors = F)

df_to_plot <-
  data.frame(samples = colnames(cpm_gene_count),
             CPM = unlist(cpm_gene_count[rownames(cpm_gene_count) %in% 
                                           "VPS45_ENST00000369130"
                                         , ]),
             genotype = df_metadata$genotype,
             stringsAsFactors = F)

### Use ENSE 1071807 that Hanwen's primer targeted region
df_to_plot <-
  data.frame(samples = colnames(cpm_gene_count),
             CPM = unlist(cpm_gene_count[rownames(cpm_gene_count) %in% 
                                           "VPS45_ENSE00001071807"
                                         , ]),
             genotype = df_metadata$genotype,
             stringsAsFactors = F)

df_to_plot <-
  data.frame(samples = colnames(cpm_gene_count),
             CPM = unlist(cpm_gene_count[rownames(cpm_gene_count) %in% 
                                           "VPS45"
                                         , ]),
             genotype = df_metadata$genotype,
             stringsAsFactors = F)

res.aov <-
  aov(CPM ~ genotype,
      data = df_to_plot)
summary(res.aov)
TukeyHSD(res.aov)

ggplot(df_to_plot,
       aes(x = genotype,
           y = CPM,
           fill = genotype)) +
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = brewer.pal(name = "Dark2", n = 3)) +
  stat_summary(aes(x = genotype,
                   y = CPM),
               fun.y = mean, geom = "point",
               size = 2, shape = 23,
               color = "red", fill = "yellow") +
  # ylim(c(0, 15)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  ggtitle("VPS45 ENST00000369130.3")

save.image(file = "make_barplot.RData")
