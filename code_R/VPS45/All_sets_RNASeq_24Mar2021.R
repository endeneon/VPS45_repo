# Siwei 25 Mar 2021
# edgeR analysis for all VPS45 RNASeq set data

# init
library(edgeR)
library(readr)

library(Rfast)
library(stringr)


# load data
gene_count_table <- read_delim("ReadsPerGene_STAR_full.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

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


# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

gene_list <- gene_count_table$Geneid
gene_count_table <- gene_count_table[, -1]
rownames(gene_count_table) <- gene_list

#### check MDS only

DGE_Raw <- DGEList(as.matrix(gene_count_table),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45", 
                                    "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                    "GG_EGFP", "GG_EGFP", "GG_EGFP")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

plotMDS(DGE_Raw)

DGE_cpm <- as.data.frame(cpm(DGE_Raw), 
                         stringsAsFactors = F)
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000285184", ]
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000136631", ]
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000118292", ]

DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000284202", ] #miR137

hist(rowmeans(as.matrix(DGE_cpm)),
     breaks = 10000)

write.table(DGE_cpm,
            file = "CPM_full_table.txt",
            row.names = T, col.names = T,
            quote = F)

## !! remove sample A11-VPS-2, A11-lnR-4, and H12-lnR-1
## !! for their abnormal locations in the MDS plot
## !! also remove 3 H6-GG samples since they cannot be used as blocking factors

gene_count_table <- gene_count_table[, -c(10, 12, 20, 
                                          26:28)]
rownames(gene_count_table) <- 
  gene_list # remember to reassign gene names after column removal

DGE_Raw <- DGEList(as.matrix(gene_count_table),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                    "AA_VPS45", "AA_VPS45", 
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", 
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = rownames(gene_count_table),
                   samples = colnames(gene_count_table),
                   remove.zeros = T)

# define samples
samples <- data.frame(genotype = c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                   "AA_VPS45", "AA_VPS45", 
                                   "AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_lncRNA", "AA_lncRNA", 
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                      cell.line = c(rep_len("A11", length.out = 11), 
                                    rep_len("H12", length.out = 11)))

design_samples <- model.matrix(~ factor(samples$cell.line) +
                                 factor(samples$genotype,
                                        levels = c("AA_EGFP", "AA_C1orf54", 
                                                   "AA_lncRNA", "AA_VPS45")))

## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_C1orf54"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_lncRNA"] >= cpm_cutoff) >= 5) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 5) 
          
          , ]

DGE_Raw <- calcNormFactors(DGE_Raw)
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

## generate different gene comparison lists

### C1orf54
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 3) # C1orf54
table(decideTests(DGE_QLF_test_results))
DGE_QLF_table <- DGE_QLF_test_results$table
DGE_QLF_table$FDR <- p.adjust(DGE_QLF_table$PValue,
                              method = "fdr")
DGE_QLF_table$Geneid <- rownames(DGE_QLF_test_results$table)
DGE_QLF_table <- merge(DGE_QLF_table,
                       gencode_v35_ENSG_Genename_final,
                       by = "Geneid")

DGE_QLF_table[DGE_QLF_table$Geneid %in% "ENSG00000118292", ]

write.table(DGE_QLF_table,
            file = "Rapid_neuron_gene_exp_table_C1orf54KD_vs_EGFP_31Mar2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")

### lncRNA
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 4) # lncRNA
table(decideTests(DGE_QLF_test_results))
DGE_QLF_table <- DGE_QLF_test_results$table
DGE_QLF_table$FDR <- p.adjust(DGE_QLF_table$PValue,
                              method = "fdr")
DGE_QLF_table$Geneid <- rownames(DGE_QLF_test_results$table)
DGE_QLF_table <- merge(DGE_QLF_table,
                       gencode_v35_ENSG_Genename_final,
                       by = "Geneid")

DGE_QLF_table[DGE_QLF_table$Gene_symbol %in% "AC244033.2", ]

write.table(DGE_QLF_table,
            file = "Rapid_neuron_gene_exp_table_lncRNAKD_vs_EGFP_31Mar2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")

### VPS45
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 5) # VPS45
table(decideTests(DGE_QLF_test_results))
DGE_QLF_table <- DGE_QLF_test_results$table
DGE_QLF_table$FDR <- p.adjust(DGE_QLF_table$PValue,
                              method = "fdr")
DGE_QLF_table$Geneid <- rownames(DGE_QLF_test_results$table)
DGE_QLF_table <- merge(DGE_QLF_table,
                       gencode_v35_ENSG_Genename_final,
                       by = "Geneid")

DGE_QLF_table[DGE_QLF_table$Gene_symbol %in% "VPS45", ]

write.table(DGE_QLF_table,
            file = "Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_31Mar2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")

#### GO analysis #####
# load data
gene_count_table <- read_delim("ReadsPerGene_STAR_full.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

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


# rename and rearrange Geneid
gene_count_table$Geneid <- str_split(gene_count_table$Geneid,
                                     pattern = "\\.",
                                     simplify = T)[, 1]
gene_count_table <- 
  gene_count_table[!duplicated(gene_count_table$Geneid), ]

gene_list <- gene_count_table$Geneid
gene_count_table <- gene_count_table[, -1]
rownames(gene_count_table) <- gene_list

gene_count_table <- gene_count_table[, -c(10, 12, 20, 
                                          26:28)]
rownames(gene_count_table) <- 
  gene_list # remember to reassign gene names after column removal

## import Entrez reference table #####
hg38_ENSG_EntrezGeneid <- read_delim("ENSG_EntrezID_hs", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

## convert ENSG to Entrez ID
hg38_ENSG_EntrezGeneid <- 
  hg38_ENSG_EntrezGeneid[!is.na(hg38_ENSG_EntrezGeneid$`NCBI gene (formerly Entrezgene) ID`), ]
hg38_ENSG_EntrezGeneid <-
  hg38_ENSG_EntrezGeneid[!duplicated(hg38_ENSG_EntrezGeneid$`NCBI gene (formerly Entrezgene) ID`), ]

Entrez_count_table <- as.data.frame(gene_count_table)

Entrez_count_table$Geneid <- rownames(Entrez_count_table)
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


# rownames(gene_count_table) <- 
#   gene_list # remember to reassign gene names after column removal

DGE_Raw <- DGEList(as.matrix(Entrez_count_table),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                    "AA_VPS45", "AA_VPS45", 
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_lncRNA", "AA_lncRNA", 
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = rownames(Entrez_count_table),
                   samples = colnames(Entrez_count_table),
                   remove.zeros = T)

# define samples
samples <- data.frame(genotype = c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                   "AA_VPS45", "AA_VPS45", 
                                   "AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_lncRNA", "AA_lncRNA", 
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                      cell.line = c(rep_len("A11", length.out = 11), 
                                    rep_len("H12", length.out = 11)))


design_samples <- model.matrix(~ factor(samples$cell.line) +
                                 factor(samples$genotype,
                                        levels = c("AA_EGFP", "AA_C1orf54", 
                                                   "AA_lncRNA", "AA_VPS45")))

## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_C1orf54"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_lncRNA"] >= cpm_cutoff) >= 5) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 5) 
          
          , ]

DGE_Raw <- calcNormFactors(DGE_Raw)
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

## generate different gene comparison lists

### C1orf54
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 3) # C1orf54
table(decideTests(DGE_QLF_test_results))

GO_results <- goana(DGE_QLF_test_results, 
                    species = "Hs",
                    FDR = 0.05,
                    trend = T,
                    plot = T)


GO_results_output <- GO_results[GO_results$Ont == "BP", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                 method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                   method = "fdr")
write.table(GO_results_output,
            file = "GO_results_BP_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "CC", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                 method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                   method = "fdr")
write.table(GO_results_output,
            file = "GO_results_CC_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "MF", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                 method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                   method = "fdr")
write.table(GO_results_output,
            file = "GO_results_MF_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

###

KEGG_results <- kegga(DGE_QLF_test_results, 
                      species = "Hs",
                      FDR = 0.05,
                      trend = T,
                      plot = T)
KEGG_results$FDR.Up <- p.adjust(KEGG_results$P.Up,
                                method = "fdr")
KEGG_results$FDR.Down <- p.adjust(KEGG_results$P.Down,
                                  method = "fdr")
write.table(KEGG_results,
            file = "KEGG_results_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")


### lncRNA
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 4) # lncRNA
table(decideTests(DGE_QLF_test_results))

GO_results <- goana(DGE_QLF_test_results, 
                    species = "Hs",
                    FDR = 0.05,
                    trend = T,
                    plot = T)


GO_results_output <- GO_results[GO_results$Ont == "BP", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_BP_lncRNAKD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "CC", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_CC_lncRNAKD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "MF", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_MF_lncRNAKD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

###

KEGG_results <- kegga(DGE_QLF_test_results, 
                      species = "Hs",
                      FDR = 0.05,
                      trend = T,
                      plot = T)
KEGG_results$FDR.Up <- p.adjust(KEGG_results$P.Up,
                                method = "fdr")
KEGG_results$FDR.Down <- p.adjust(KEGG_results$P.Down,
                                  method = "fdr")
write.table(KEGG_results,
            file = "KEGG_results_lncRNAKD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")


### VPS45
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 5) # VPS45
table(decideTests(DGE_QLF_test_results))

GO_results <- goana(DGE_QLF_test_results, 
                    species = "Hs",
                    FDR = 0.05,
                    trend = T,
                    plot = T)


GO_results_output <- GO_results[GO_results$Ont == "BP", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_BP_VPS45KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "CC", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_CC_VPS45KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "MF", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_MF_VPS45KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

###

KEGG_results <- kegga(DGE_QLF_test_results, 
                      species = "Hs",
                      FDR = 0.05,
                      trend = T,
                      plot = T)
KEGG_results$FDR.Up <- p.adjust(KEGG_results$P.Up,
                                method = "fdr")
KEGG_results$FDR.Down <- p.adjust(KEGG_results$P.Down,
                                  method = "fdr")
write.table(KEGG_results,
            file = "KEGG_results_VPS45KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")


######

#AAvsGG

Entrez_count_table <- Entrez_count_table[, -c(1:3, 7:16, 20:25)]

rownames(Entrez_count_table) <- 
  gene_list # remember to reassign gene names after column removal

DGE_Raw <- DGEList(as.matrix(Entrez_count_table),
                   group = factor(c(#"AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    # "AA_lncRNA", "AA_lncRNA", "AA_lncRNA", 
                                    # "AA_VPS45", "AA_VPS45", 
                                    # "AA_C1orf54", "AA_C1orf54", "AA_C1orf54", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    # "AA_VPS45", "AA_VPS45", 
                                    # "AA_lncRNA", "AA_lncRNA", "AA_lncRNA",
                                    "GG_EGFP", "GG_EGFP", "GG_EGFP")),
                   genes = rownames(Entrez_count_table),
                   samples = colnames(Entrez_count_table),
                   remove.zeros = T)

# define samples
samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "GG_EGFP", "GG_EGFP", "GG_EGFP"))

design_samples <- model.matrix(~ factor(samples$genotype,
                                        levels = c("AA_EGFP", "GG_EGFP")))

## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "GG_EGFP"] >= cpm_cutoff) >= 3)
          
          , ]

DGE_Raw <- calcNormFactors(DGE_Raw)
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

## generate different gene comparison lists

### GG vs AA
DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 2) # GG vs AA
table(decideTests(DGE_QLF_test_results))

GO_results <- goana(DGE_QLF_test_results, 
                    species = "Hs",
                    FDR = 0.05,
                    trend = T,
                    plot = T)


GO_results_output <- GO_results[GO_results$Ont == "BP", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_BP_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "CC", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_CC_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

GO_results_output <- GO_results[GO_results$Ont == "MF", ]
GO_results_output$FDR.Up <- p.adjust(GO_results_output$P.Up,
                                     method = "fdr")
GO_results_output$FDR.Down <- p.adjust(GO_results_output$P.Down,
                                       method = "fdr")
write.table(GO_results_output,
            file = "GO_results_MF_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

###

KEGG_results <- kegga(DGE_QLF_test_results, 
                      species = "Hs",
                      FDR = 0.05,
                      trend = T,
                      plot = T)
KEGG_results$FDR.Up <- p.adjust(KEGG_results$P.Up,
                                method = "fdr")
KEGG_results$FDR.Down <- p.adjust(KEGG_results$P.Down,
                                  method = "fdr")
write.table(KEGG_results,
            file = "KEGG_results_C1orf54KD_vs_EGFP.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")
