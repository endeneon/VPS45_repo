# Siwei 01 Apr 2021
# DeSeq2 analysis for all VPS45 RNASeq set data

# init
library(DESeq2)
library(BiocParallel)

library(readr)

library(Rfast)
library(stringr)

# init multicore
register(MulticoreParam(8))

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
rownames(samples) <- colnames(gene_count_table) # make consistent

dds <- DESeqDataSetFromMatrix(countData = gene_count_table,
                              colData = samples,
                              design = ~ cell.line + genotype)

ddsMF <- dds[rowSums(counts(dds)) > 2, ]
ddsMF <- DESeq(ddsMF)
resMFgenotype <- results(ddsMF,
                         contrast = c("genotype", "AA_VPS45", "AA_EGFP"))

summary(resMFgenotype)

resMFgenotype_output <- as.data.frame(resMFgenotype)
resMFgenotype_output <- resMFgenotype_output[!is.na(resMFgenotype_output$padj), ]
sum(resMFgenotype_output$padj < 0.05)
resMFgenotype_output$Geneid <- rownames(resMFgenotype_output)

resMFgenotype_output <- merge(resMFgenotype_output,
                              gencode_v35_ENSG_Genename_final,
                              by = "Geneid")
write.table(resMFgenotype_output,
            file = "DeSeq2_Rapid_neuron_gene_exp_table_VPS45KD_vs_EGFP_01Apr2021.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")
