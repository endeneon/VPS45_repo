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


## !! remove sample A11-VPS-2, A11-lnR-4, and H12-lnR-1
## !! for their abnormal locations in the MDS plot
## !! also remove 3 H6-GG samples since they cannot be used as blocking factors

gene_count_table <- gene_count_table[, -c(10, 12, 20, 
                                          26:28)]
rownames(gene_count_table) <- 
  gene_list # remember to reassign gene names after column removal

## analyze separate data sets
## VPS45

gene_count_table_sub <- gene_count_table[, c(4:6, 10, 11, 
                                             15:17, 20:22)]
rownames(gene_count_table_sub) <- gene_list


DGE_Raw <- DGEList(as.matrix(gene_count_table_sub),
                   group = factor(c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", 
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_VPS45", "AA_VPS45", "AA_VPS45")),
                   genes = rownames(gene_count_table_sub),
                   samples = colnames(gene_count_table_sub),
                   remove.zeros = T)

# define samples
samples <- data.frame(genotype = c("AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", 
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_VPS45", "AA_VPS45", "AA_VPS45"),
                      cell.line = c(rep_len("A11", length.out = 5), 
                                    rep_len("H12", length.out = 6)))

design_samples <- model.matrix(~ factor(samples$cell.line) +
                                 factor(samples$genotype,
                                        levels = c("AA_EGFP", "AA_VPS45")))

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

DGE_QLF_test_results <- glmQLFTest(DGE_QLF,
                                   coef = 3) # VPS45
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


## C1orf54

gene_count_table_sub <- gene_count_table[, c(1:6, 12:17)]
rownames(gene_count_table_sub) <- gene_list

DGE_Raw <- DGEList(as.matrix(gene_count_table_sub),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP")),
                   genes = rownames(gene_count_table_sub),
                   samples = colnames(gene_count_table_sub),
                   remove.zeros = T)

# define samples
samples <- data.frame(genotype = c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP"),
                      cell.line = c(rep_len("A11", length.out = 6), 
                                    rep_len("H12", length.out = 6)))

design_samples <- model.matrix(~ factor(samples$cell.line) +
                                 factor(samples$genotype,
                                        levels = c("AA_EGFP", "AA_C1orf54")))

## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
cpm_gene_count <- as.data.frame(cpm_gene_count,
                                stringsAsFactors = F)
cpm_gene_count$Geneid <- DGE_Raw$genes$genes
cpm_gene_count[cpm_gene_count$Geneid %in% "ENSG00000118292", ]
cpm_gene_count[cpm_gene_count$Geneid %in% "ENSG00000118292", ]


## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 6) 
          
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



## lncRNA
gene_count_table_sub <- gene_count_table[, c(4:9, 15:19)]
rownames(gene_count_table_sub) <- gene_list

DGE_Raw <- DGEList(as.matrix(gene_count_table_sub),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP")),
                   genes = rownames(gene_count_table_sub),
                   samples = colnames(gene_count_table_sub),
                   remove.zeros = T)

# define samples
samples <- data.frame(genotype = c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                   "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                   "AA_EGFP", "AA_EGFP", "AA_EGFP"),
                      cell.line = c(rep_len("A11", length.out = 6), 
                                    rep_len("H12", length.out = 6)))

design_samples <- model.matrix(~ factor(samples$cell.line) +
                                 factor(samples$genotype,
                                        levels = c("AA_EGFP", "AA_C1orf54")))

## check cpm histogram and filter values
cpm_gene_count <- cpm(DGE_Raw)
hist(cpm(DGE_Raw), 
     breaks = 100000,
     xlim = c(0, 1))
cpm_gene_count <- as.data.frame(cpm_gene_count,
                                stringsAsFactors = F)
cpm_gene_count$Geneid <- DGE_Raw$genes$genes
cpm_gene_count[cpm_gene_count$Geneid %in% "ENSG00000118292", ]
cpm_gene_count[cpm_gene_count$Geneid %in% "ENSG00000118292", ]


## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_VPS45"] >= cpm_cutoff) >= 6) 
          
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
