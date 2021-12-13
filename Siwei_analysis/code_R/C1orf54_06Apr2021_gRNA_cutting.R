# Siwei 06 Apr 2021
# edgeR analysis for Alena's C1orf54 set
# extract cpm and DE results including the C1orf54 exons 
# before the gRNA cutting site

# init
library(edgeR)
library(readr)

library(Rfast)
library(stringr)

library(gplots)

# load data
gene_count_table <- read_delim("ReadsPerGene_STAR_full.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

C1orf54_gRNA_input <- read_delim("C1orf54_gRNA_input.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
C1orf54_gene_list <- C1orf54_gRNA_input$Geneid
C1orf54_gRNA_input <- C1orf54_gRNA_input[-c(1:6)]
rownames(C1orf54_gRNA_input) <- C1orf54_gene_list

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
# cpm_gene_count[cpm_gene_count$Geneid %in% "ENSG00000118292", ]
# cpm_gene_count[cpm_gene_count$Geneid %in% "ENSG00000118292", ]


## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 1

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$genotype %in% "AA_EGFP"] >= cpm_cutoff) >= 6) |
            (rowSums(cpm_gene_count[, samples$genotype %in% "AA_C1orf54"] >= cpm_cutoff) >= 6) 
          
          , ]

## ! obtain new gene table after cpm cut-off
gene_count_table_post_cpm_cutoff <-
  gene_count_table_sub[rownames(gene_count_table_sub) %in% DGE_Raw$genes$genes, ]
# add two C1orf54 fragment back
# extract gene ids from gene_count_table_post_cpm_cutoff and save first
gene_count_table_post_cpm_cutoff_id <- DGE_Raw$genes$genes
gene_count_table_post_cpm_cutoff <-
  data.frame(rbind(gene_count_table_post_cpm_cutoff,
                   C1orf54_gRNA_input),
             stringsAsFactors = F, 
             row.names = c(gene_count_table_post_cpm_cutoff_id,
                           rownames(C1orf54_gRNA_input)))

## We have re-assembled gene count table now, 
## we can re-run EdgeR again
## note we should not apply CPM cut-off this time to avoid
## accidentally removal of the two C1orf54 fragments

DGE_Raw <- DGEList(as.matrix(gene_count_table_post_cpm_cutoff),
                   group = factor(c("AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP",
                                    "AA_C1orf54", "AA_C1orf54", "AA_C1orf54",
                                    "AA_EGFP", "AA_EGFP", "AA_EGFP")),
                   genes = rownames(gene_count_table_post_cpm_cutoff),
                   samples = colnames(gene_count_table_post_cpm_cutoff),
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

DGE_Raw <- calcNormFactors(DGE_Raw,
                           method = "TMM")
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

## ! note here we need to put all.x = T to avoid removal of the two C1orf54 frags
DGE_QLF_table <- merge(DGE_QLF_table,
                       gencode_v35_ENSG_Genename_final,
                       all.x = T,
                       by = "Geneid")

DGE_QLF_table[DGE_QLF_table$Geneid %in% "ENSG00000118292", ]

write.table(DGE_QLF_table,
            file = "Gene_exp_table_C1orf54KD_vs_EGFP_inc_C1orf54_frags_06Apr2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")


cpm_gene_count <- merge(cpm_gene_count,
                        gencode_v35_ENSG_Genename_final,
                        all.x = T,
                        by = "Geneid")
write.table(cpm_gene_count,
            file = "CPM_table_DE_genes_only_C1orf54KD_vs_EGFP_inc_C1orf54_frags_06Apr2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")
