# calculate the correlation of A11 vs H12 lines

# init
library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)

library(ggplot2)
library(gplots)
library(RColorBrewer)

# load data
gene_count_table <- read_delim("ReadsPerGene_STAR.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
gene_count_table$...17 <- NULL

# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

count_group_1 <- gene_count_table[, 2:7]
rownames(count_group_1) <- gene_count_table$Geneid


## make the raw DGEList from the count matrix
DGE_Raw_1 <- DGEList(as.matrix(count_group_1),
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = gene_count_table$Geneid,
                   samples = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                               "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                   remove.zeros = T)

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count_1 <- cpm(DGE_Raw_1)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw_1 <- 
  DGE_Raw_1[(rowSums(cpm_gene_count_1[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 2) |
            (rowSums(cpm_gene_count_1[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 2) 
          , ]

DGE_Raw_1 <- calcNormFactors(DGE_Raw_1,
                             method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw_1 <- estimateDisp(DGE_Raw_1, 
                          design = design_samples,
                          robust = T)


DGE_QLF_1 <- glmQLFit(DGE_Raw_1,
                      design = design_samples,
                      robust = T)
DGE_QLF_1 <- glmQLFTest(DGE_QLF_1,
                        coef = 2) 

####### calc group 2

count_group_2 <- gene_count_table[, 8:13]
rownames(count_group_2) <- gene_count_table$Geneid

## make the raw DGEList from the count matrix
DGE_Raw_2 <- DGEList(as.matrix(count_group_2),
                     group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                      "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                     genes = gene_count_table$Geneid,
                     samples = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                 "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                     remove.zeros = T)

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count_2 <- cpm(DGE_Raw_2)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw_2 <- 
  DGE_Raw_2[(rowSums(cpm_gene_count_2[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 2) |
              (rowSums(cpm_gene_count_2[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 2) 
            , ]

DGE_Raw_2 <- calcNormFactors(DGE_Raw_2,
                             method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw_2 <- estimateDisp(DGE_Raw_2, 
                          design = design_samples,
                          robust = T)

DGE_QLF_2 <- glmQLFit(DGE_Raw_2,
                      design = design_samples,
                      robust = T)
DGE_QLF_2 <- glmQLFTest(DGE_QLF_2,
                        coef = 2) 

##### calc correlation
calc_table_g1 <- DGE_QLF_1$table
calc_table_g1$FDR <- p.adjust(calc_table_g1$PValue,
                              method = "fdr")
# calc_table_g1 <- calc_table_g1[calc_table_g1$FDR < 0.05, ]
# calc_table_g1 <- calc_table_g1[calc_table_g1$PValue < 0.05, ]
calc_table_g1$Geneid <- rownames(calc_table_g1)
# calc_table_g1 <- calc_table_g1[order(calc_table_g1$PValue), ]
# calc_table_g1 <- calc_table_g1[1:5000, ]

calc_table_g2 <- DGE_QLF_2$table
calc_table_g2$FDR <- p.adjust(calc_table_g2$PValue,
                              method = "fdr")
# calc_table_g2 <- calc_table_g1[calc_table_g2$FDR < 0.05, ]
# calc_table_g2 <- calc_table_g1[calc_table_g2$PValue < 0.05, ]
calc_table_g2$Geneid <- rownames(calc_table_g2)
# calc_table_g2 <- calc_table_g2[order(calc_table_g2$PValue), ]
# calc_table_g2 <- calc_table_g2[1:5000, ]

###
# calc_table_merged <-
#   merge(x = calc_table_g1,
#         y = calc_table_g2,
#         by = "Geneid")

calc_table_common_index <-
  c(calc_table_g1$Geneid[calc_table_g1$FDR < 0.05],
    calc_table_g2$Geneid[calc_table_g2$FDR < 0.05])

calc_table_common_index <-
  calc_table_common_index[!duplicated(calc_table_common_index)]
calc_table_common_index <- 
  calc_table_common_index[calc_table_common_index %in%
                            calc_table_g1$Geneid]
calc_table_common_index <- 
  calc_table_common_index[calc_table_common_index %in%
                            calc_table_g2$Geneid]

calc_table_merged <-
  data.frame(Geneid = calc_table_common_index,
             logFC.x = calc_table_g1$logFC[calc_table_g1$Geneid %in% 
                                             calc_table_common_index],
             logFC.y = calc_table_g2$logFC[calc_table_g2$Geneid %in% 
                                             calc_table_common_index],
             stringsAsFactors = F)


cor.test(x = calc_table_merged$logFC.x,
         y = calc_table_merged$logFC.y,
         method = "pearson",
         alternative = "t")

cor.test(x = calc_table_merged$logFC.x,
         y = calc_table_merged$logFC.y,
         method = "spearman",
         alternative = "t")

ggplot(calc_table_merged,
       aes(x = logFC.x,
           y = logFC.y)) +
  geom_point(size = 0.3) +
  stat_smooth(method = "glm") +
  xlab("logFC.AA_EGFP") +
  ylab("logFC.AA_VPS45") +
  theme_classic()

##### C1orf54
# load data
gene_count_table <- read_delim("ReadsPerGene_STAR_RNASeq_C1orf54_VPS45_lncRNA_GG.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
# gene_count_table$...17 <- NULL

# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

count_group_1 <- gene_count_table[, c(5:7, 2:4)]
rownames(count_group_1) <- gene_count_table$Geneid


## make the raw DGEList from the count matrix
DGE_Raw_1 <- DGEList(as.matrix(count_group_1),
                     group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                      "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                     genes = gene_count_table$Geneid,
                     samples = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                 "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                     remove.zeros = T)

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count_1 <- cpm(DGE_Raw_1)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw_1 <- 
  DGE_Raw_1[(rowSums(cpm_gene_count_1[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 2) |
              (rowSums(cpm_gene_count_1[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 2) 
            , ]

DGE_Raw_1 <- calcNormFactors(DGE_Raw_1,
                             method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw_1 <- estimateDisp(DGE_Raw_1, 
                          design = design_samples,
                          robust = T)


DGE_QLF_1 <- glmQLFit(DGE_Raw_1,
                      design = design_samples,
                      robust = T)
DGE_QLF_1 <- glmQLFTest(DGE_QLF_1,
                        coef = 2) 

####### calc group 2

count_group_2 <- gene_count_table[, c(17:19, 14:16)]
rownames(count_group_2) <- gene_count_table$Geneid

## make the raw DGEList from the count matrix
DGE_Raw_2 <- DGEList(as.matrix(count_group_2),
                     group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                      "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                     genes = gene_count_table$Geneid,
                     samples = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                 "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                     remove.zeros = T)

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count_2 <- cpm(DGE_Raw_2)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw_2 <- 
  DGE_Raw_2[(rowSums(cpm_gene_count_2[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 2) |
              (rowSums(cpm_gene_count_2[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 2) 
            , ]

DGE_Raw_2 <- calcNormFactors(DGE_Raw_2,
                             method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw_2 <- estimateDisp(DGE_Raw_2, 
                          design = design_samples,
                          robust = T)

DGE_QLF_2 <- glmQLFit(DGE_Raw_2,
                      design = design_samples,
                      robust = T)
DGE_QLF_2 <- glmQLFTest(DGE_QLF_2,
                        coef = 2) 

##### calc correlation
calc_table_g1 <- DGE_QLF_1$table
calc_table_g1$FDR <- p.adjust(calc_table_g1$PValue,
                              method = "fdr")
# calc_table_g1 <- calc_table_g1[calc_table_g1$FDR < 0.05, ]
# calc_table_g1 <- calc_table_g1[calc_table_g1$PValue < 0.05, ]
calc_table_g1$Geneid <- rownames(calc_table_g1)
# calc_table_g1 <- calc_table_g1[order(calc_table_g1$PValue), ]
# calc_table_g1 <- calc_table_g1[1:5000, ]

calc_table_g2 <- DGE_QLF_2$table
calc_table_g2$FDR <- p.adjust(calc_table_g2$PValue,
                              method = "fdr")
# calc_table_g2 <- calc_table_g1[calc_table_g2$FDR < 0.05, ]
# calc_table_g2 <- calc_table_g1[calc_table_g2$PValue < 0.05, ]
calc_table_g2$Geneid <- rownames(calc_table_g2)
# calc_table_g2 <- calc_table_g2[order(calc_table_g2$PValue), ]
# calc_table_g2 <- calc_table_g2[1:5000, ]

###
# calc_table_merged <-
#   merge(x = calc_table_g1,
#         y = calc_table_g2,
#         by = "Geneid")

calc_table_common_index <-
  c(calc_table_g1$Geneid[calc_table_g1$FDR < 0.05],
    calc_table_g2$Geneid[calc_table_g2$FDR < 0.05])

calc_table_common_index <-
  calc_table_common_index[!duplicated(calc_table_common_index)]
calc_table_common_index <- 
  calc_table_common_index[calc_table_common_index %in%
                            calc_table_g1$Geneid]
calc_table_common_index <- 
  calc_table_common_index[calc_table_common_index %in%
                            calc_table_g2$Geneid]

calc_table_merged <-
  data.frame(Geneid = calc_table_common_index,
             logFC.x = calc_table_g1$logFC[calc_table_g1$Geneid %in% 
                                             calc_table_common_index],
             logFC.y = calc_table_g2$logFC[calc_table_g2$Geneid %in% 
                                             calc_table_common_index],
             stringsAsFactors = F)


cor.test(x = calc_table_merged$logFC.x,
         y = calc_table_merged$logFC.y,
         method = "pearson",
         alternative = "t")

# cor.test(x = calc_table_merged$logFC.x,
#          y = calc_table_merged$logFC.y,
#          method = "spearman",
#          alternative = "t")

ggplot(calc_table_merged,
       aes(x = logFC.x,
           y = logFC.y)) +
  geom_point(size = 0.3) +
  stat_smooth(method = "glm") +
  xlab("logFC.AA_EGFP") +
  ylab("logFC.AA_C1orf54") +
  theme_classic()


####################lncRNA
# load data
gene_count_table <- read_delim("ReadsPerGene_STAR_RNASeq_C1orf54_VPS45_lncRNA_GG.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
# gene_count_table$...17 <- NULL

# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

count_group_1 <- gene_count_table[, c(5:7, 8:10)]
rownames(count_group_1) <- gene_count_table$Geneid


## make the raw DGEList from the count matrix
DGE_Raw_1 <- DGEList(as.matrix(count_group_1),
                     group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                      "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                     genes = gene_count_table$Geneid,
                     samples = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                 "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                     remove.zeros = T)

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count_1 <- cpm(DGE_Raw_1)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw_1 <- 
  DGE_Raw_1[(rowSums(cpm_gene_count_1[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 2) |
              (rowSums(cpm_gene_count_1[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 2) 
            , ]

DGE_Raw_1 <- calcNormFactors(DGE_Raw_1,
                             method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw_1 <- estimateDisp(DGE_Raw_1, 
                          design = design_samples,
                          robust = T)


DGE_QLF_1 <- glmQLFit(DGE_Raw_1,
                      design = design_samples,
                      robust = T)
DGE_QLF_1 <- glmQLFTest(DGE_QLF_1,
                        coef = 2) 

####### calc group 2

count_group_2 <- gene_count_table[, c(17:19, 18:20)]
rownames(count_group_2) <- gene_count_table$Geneid

## make the raw DGEList from the count matrix
DGE_Raw_2 <- DGEList(as.matrix(count_group_2),
                     group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                      "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                     genes = gene_count_table$Geneid,
                     samples = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                 "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                     remove.zeros = T)

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count_2 <- cpm(DGE_Raw_2)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw_2 <- 
  DGE_Raw_2[(rowSums(cpm_gene_count_2[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 2) |
              (rowSums(cpm_gene_count_2[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 2) 
            , ]

DGE_Raw_2 <- calcNormFactors(DGE_Raw_2,
                             method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw_2 <- estimateDisp(DGE_Raw_2, 
                          design = design_samples,
                          robust = T)

DGE_QLF_2 <- glmQLFit(DGE_Raw_2,
                      design = design_samples,
                      robust = T)
DGE_QLF_2 <- glmQLFTest(DGE_QLF_2,
                        coef = 2) 

##### calc correlation
calc_table_g1 <- DGE_QLF_1$table
calc_table_g1$FDR <- p.adjust(calc_table_g1$PValue,
                              method = "fdr")
# calc_table_g1 <- calc_table_g1[calc_table_g1$FDR < 0.05, ]
# calc_table_g1 <- calc_table_g1[calc_table_g1$PValue < 0.05, ]
calc_table_g1$Geneid <- rownames(calc_table_g1)
# calc_table_g1 <- calc_table_g1[order(calc_table_g1$PValue), ]
# calc_table_g1 <- calc_table_g1[1:5000, ]

calc_table_g2 <- DGE_QLF_2$table
calc_table_g2$FDR <- p.adjust(calc_table_g2$PValue,
                              method = "fdr")
# calc_table_g2 <- calc_table_g1[calc_table_g2$FDR < 0.05, ]
# calc_table_g2 <- calc_table_g1[calc_table_g2$PValue < 0.05, ]
calc_table_g2$Geneid <- rownames(calc_table_g2)
# calc_table_g2 <- calc_table_g2[order(calc_table_g2$PValue), ]
# calc_table_g2 <- calc_table_g2[1:5000, ]

###
# calc_table_merged <-
#   merge(x = calc_table_g1,
#         y = calc_table_g2,
#         by = "Geneid")



calc_table_common_index <-
  c(calc_table_g1$Geneid[calc_table_g1$FDR < 0.05],
    calc_table_g2$Geneid[calc_table_g2$FDR < 0.05])

calc_table_common_index <-
  calc_table_common_index[!duplicated(calc_table_common_index)]
calc_table_common_index <- 
  calc_table_common_index[calc_table_common_index %in%
                            calc_table_g1$Geneid]
calc_table_common_index <- 
  calc_table_common_index[calc_table_common_index %in%
                            calc_table_g2$Geneid]

calc_table_merged <-
  data.frame(Geneid = calc_table_common_index,
             logFC.x = calc_table_g1$logFC[calc_table_g1$Geneid %in% 
                                             calc_table_common_index],
             logFC.y = calc_table_g2$logFC[calc_table_g2$Geneid %in% 
                                             calc_table_common_index],
             stringsAsFactors = F)

calc_table_common_index <- merge(x = calc_table_g1,
                                 y = calc_table_g2,
                                 by = "Geneid")
calc_table_common_index <- 
  calc_table_common_index[sample(x = nrow(calc_table_common_index),
                                 size = 2500), 
                          ]

cor.test(x = calc_table_common_index$logFC.x,
         y = calc_table_common_index$logFC.y,
         method = "pearson",
         alternative = "t")

# cor.test(x = calc_table_merged$logFC.x,
#          y = calc_table_merged$logFC.y,
#          method = "spearman",
#          alternative = "t")

### calc random merged samples




ggplot(calc_table_merged,
       aes(x = logFC.x,
           y = logFC.y)) +
  geom_point(size = 0.3) +
  stat_smooth(method = "glm") +
  xlab("logFC.AA_EGFP") +
  ylab("logFC.AA_lncRNA") +
  theme_classic()

