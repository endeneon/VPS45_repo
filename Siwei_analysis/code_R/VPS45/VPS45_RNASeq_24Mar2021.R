# Siwei 24 Mar 2021
# edgeR analysis for VPS45 RNASeq data

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
gene_count_table$X17 <- NULL

# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

#
## import Entrez reference table #####
hg38_ENSG_EntrezGeneid <- read_delim("ENSG_EntrezID_hs", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

## convert ENSG to Entrez ID
hg38_ENSG_EntrezGeneid <- 
  hg38_ENSG_EntrezGeneid[!is.na(hg38_ENSG_EntrezGeneid$`NCBI gene (formerly Entrezgene) ID`), ]
hg38_ENSG_EntrezGeneid <-
  hg38_ENSG_EntrezGeneid[!duplicated(hg38_ENSG_EntrezGeneid$`NCBI gene (formerly Entrezgene) ID`), ]

Entrez_count_table <- as.data.frame(gene_count_table)

Entrez_count_table <- merge(Entrez_count_table,
                            hg38_ENSG_EntrezGeneid,
                            by.x = "Geneid",
                            by.y = colnames(hg38_ENSG_EntrezGeneid)[1])
rownames(Entrez_count_table) <- Entrez_count_table$`NCBI gene (formerly Entrezgene) ID`
gene_list <- rownames(Entrez_count_table)
Entrez_count_table$Geneid <- NULL
Entrez_count_table$`Gene name` <- NULL
Entrez_count_table$`NCBI gene (formerly Entrezgene) ID` <- NULL
Entrez_count_table$`Chromosome/scaffold name` <- NULL
Entrez_count_table$`Gene start (bp)` <- NULL
Entrez_count_table$`Gene end (bp)` <- NULL
Entrez_count_table$Strand <- NULL

## remove GG columns since GG columns are from a separate line and is not
## directly estimable for coefficient.
## will calculate it later
Entrez_count_table <- Entrez_count_table[ , -c(13:15)]

## make the raw DGEList from the count matrix
DGE_Raw <- DGEList(as.matrix(Entrez_count_table),
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = rownames(Entrez_count_table),
                   samples = colnames(Entrez_count_table),
                   remove.zeros = T)

DGE_cpm <- cpm(DGE_Raw)
DGE_cpm <- as.data.frame(DGE_cpm, 
                         stringsAsFactors = F)
# DGE_cpm[rownames(DGE_cpm) %in% "11311", ]

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                      cell.line = c("A11", "A11", "A11", "A11", "A11", "A11", 
                                    "H12", "H12", "H12", "H12", "H12", "H12"))

design_samples <- model.matrix(~ #0 + # note the "0" here
                                 factor(samples$cell.line) +
                                 factor(samples$genotype))

cpm_gene_count <- cpm(DGE_Raw)
hist(cpm_gene_count, 
     breaks = 100000,
     xlim = c(0, 1))
## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 3) 
          , ]

DGE_Raw <- calcNormFactors(DGE_Raw,
                           method = "none") #c("TMM","TMMwsp","RLE","upperquartile","none")
DGE_Raw <- estimateDisp(DGE_Raw, 
                        design = design_samples,
                        robust = T)

## make MDS and BCV plots to observe the clustering pattern of samples
## with the cell-type element regressed out as blocking factor.
plotMDS(DGE_Raw)
plotBCV(DGE_Raw)

DGE_QLF <- glmQLFit(DGE_Raw,
                    design = design_samples,
                    robust = T)
DGE_QLF <- glmQLFTest(DGE_QLF,
                      coef = 3) # refer to the structure of the design matrix
# DGE_QLF <- glmQLFTest(DGE_QLF,
#                       contrast = c(-1, 0, 1)) # refer to the structure of the design matrix
table(decideTests(DGE_QLF))

GO_results <- goana(DGE_QLF, 
                    species = "Hs",
                    FDR = 0.05,
                    trend = T,
                    plot = T)


GO_results_BP <- GO_results[GO_results$Ont == "BP", ]
GO_results_BP$FDR.Up <- p.adjust(GO_results_BP$P.Up,
                                 method = "fdr")
GO_results_BP$FDR.Down <- p.adjust(GO_results_BP$P.Down,
                                   method = "fdr")
write.table(GO_results_BP,
            file = "GO_results_BP_VPS45KD_vs_EGFP_2_lines.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")


##### Not using edgeR built-in GO analysis #####
## generate gene expression list only

gene_count_table <- read_delim("ReadsPerGene_STAR.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
gene_count_table$X17 <- NULL
# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

gene_list <- gene_count_table$Geneid
gene_count_table$Geneid <- NULL
gene_count_table <- gene_count_table[, -c(13:15)] ## remove GG columns 
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
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

DGE_cpm <- cpm(DGE_Raw)
DGE_cpm <- as.data.frame(DGE_cpm, 
                         stringsAsFactors = F)
# DGE_cpm[rownames(DGE_cpm) %in% "11311", ]

samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                      cell.line = c("A11", "A11", "A11", "A11", "A11", "A11", 
                                    "H12", "H12", "H12", "H12", "H12", "H12"))

# design_samples <- model.matrix(~ factor(samples$cell.line) +
#                                  factor(samples$genotype))

design_samples <- model.matrix(~ factor(samples$genotype))

cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 3) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 3) 
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

DGE_QLF <- glmQLFit(DGE_Raw,
                    design = design_samples,
                    robust = T)
DGE_QLF <- glmQLFTest(DGE_QLF,
                      coef = 2) # refer to the structure of the design matrix
# DGE_QLF <- glmQLFTest(DGE_QLF,
#                       contrast = c(-1, 0, 1)) # refer to the structure of the design matrix
table(decideTests(DGE_QLF))

####
## method = "none"
# -1     0     1 
# 2174 12648  2974 
## method = "TMM" (default)
# -1     0     1 
# 2502 12973  2321 
## method = "RLE"
# -1     0     1 
# 2520 12981  2295 
## method = "upperquartile"
# -1     0     1 
# 2334 12895  2567 
####

DGE_QLF_table <- DGE_QLF$table
DGE_QLF_table$FDR <- p.adjust(DGE_QLF_table$PValue,
                              method = "fdr")
DGE_QLF_table$Geneid <- rownames(DGE_QLF_table)

DGE_QLF_table <- merge(DGE_QLF_table,
                       gencode_v35_ENSG_Genename_final,
                       by = "Geneid")

write.table(DGE_QLF_table,
            file = "Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_25Mar2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")

### misc
DGE_cpm_merged <- merge(DGE_cpm, gencode_v35_ENSG_Genename_final,
                        by = "Geneid")
DGE_cpm_merged[DGE_cpm_merged$Gene_symbol %in% "C1orf54", ]
DGE_cpm_merged[DGE_cpm_merged$Gene_symbol %in% "AC244033.2", ]
DGE_cpm_merged[DGE_cpm_merged$Gene_symbol %in% "VPS45", ]


###

###### for test MDS plot only #####
gene_count_table <- read_delim("ReadsPerGene_STAR.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
gene_count_table$X17 <- NULL

# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]


rownames(gene_count_table) <- gene_count_table$Geneid
gene_list <- gene_count_table$Geneid
gene_count_table <- gene_count_table[, -1]
rownames(gene_count_table) <- gene_list


DGE_Raw <- DGEList(as.matrix(gene_count_table),
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "GG_EGFP", "GG_EGFP", "GG_EGFP")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

plotMDS(DGE_Raw)

DGE_cpm <- cpm(DGE_Raw)




gene_count_table <- read_delim("featurecounts_gencode_v30_by_gene_count.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE, 
                               skip = 1)
gene_count_table <- gene_count_table[, -(2:6)]
colnames(gene_count_table) <- 
  str_replace_all(colnames(gene_count_table),
                  pattern = "Aligned\\.sortedByCoord\\.out\\.bam",
                  replacement = "")

gene_count_table$X17 <- NULL
gene_list <- gene_count_table$Geneid
gene_count_table$Geneid <- NULL
rownames(gene_count_table) <- gene_list

DGE_Raw <- DGEList(as.matrix(gene_count_table),
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "GG_EGFP", "GG_EGFP", "GG_EGFP")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

plotMDS(DGE_Raw)

DGE_cpm <- cpm(DGE_Raw)
DGE_cpm <- as.data.frame(DGE_cpm,
                         stringsAsFactors = F)
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000285184.2", ]
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000118292.9", ]


DGE_cpm_merged <- merge(DGE_cpm, gencode_v35_ENSG_Genename_final,
                        by = "Geneid")
DGE_cpm_merged[DGE_cpm_merged$Gene_symbol %in% "C1orf54", ]
DGE_cpm_merged[DGE_cpm_merged$Gene_symbol %in% "AC244033.2", ]
DGE_cpm_merged[DGE_cpm_merged$Gene_symbol %in% "VPS45", ]

gencode_v35_ENSG_Genename_final[gencode_v35_ENSG_Genename_final$Gene_symbol %in% "C1orf54", ]


###### plot PCA use fviz#####
gene_count_table <- read_delim("ReadsPerGene_STAR.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
gene_count_table$X17 <- NULL

# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]


rownames(gene_count_table) <- gene_count_table$Geneid
gene_list <- gene_count_table$Geneid
gene_count_table <- gene_count_table[, -1]
rownames(gene_count_table) <- gene_list


DGE_Raw <- DGEList(as.matrix(gene_count_table),
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "GG_EGFP", "GG_EGFP", "GG_EGFP")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

plotMDS(DGE_Raw)

DGE_cpm <- cpm(DGE_Raw)
