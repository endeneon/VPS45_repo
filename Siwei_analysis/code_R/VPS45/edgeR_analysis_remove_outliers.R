# Siwei 16 Oct 2020
# edgeR analysis for Zhiping's RNASeq data
# reanalyse by removing N40 V1+V2 as outliers

# init
library(edgeR)
library(readr)
library(readxl)
library(Rfast)
library(ggplot2)
library(stringr)

##
# load data
## read in gene count matrix
gene_count_table <- readRDS(file = "hg38_gene_counts.rds")
## import metadata
samples <- read_excel("samples.xlsx")
## import Entrez reference table
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

## remove outlier data
# gene_count_table <- gene_count_table[, -c(10:11)]
Entrez_count_table <- Entrez_count_table[, -c(10:11)]
samples <- samples[-c(10:11), ]

## make the raw DGEList from the count matrix
DGE_Raw <- DGEList(as.matrix(Entrez_count_table),
                   group = as.factor(samples$treatment),
                   genes = rownames(Entrez_count_table),
                   samples = colnames(Entrez_count_table),
                   remove.zeros = T)

# ## make the raw DGEList from the count matrix
# DGE_Raw <- DGEList(as.matrix(gene_count_table),
#                    group = as.factor(samples$treatment),
#                    genes = rownames(gene_count_table),
#                    samples = colnames(gene_count_table),
#                    remove.zeros = T)

## make design matrix. Since (assumingly) we want to regress out the cell
## type as the blocking factor, I placed the cell type at the last of the 
## list. If you want to compare the cell type effect whilst regressing out
# the treatment as the blocking factor you need to switch the position of
## 'samples$treatment' and 'samples$cell'
design_samples <- model.matrix(~ 0 + # note the "0" here
                                 factor(samples$treatment) +
                                 factor(samples$cell))

## check the structure of the design matrix
head(design_samples)

## calculate cpm for optional filter-out of low-expression genes
cpm_gene_count <- cpm(DGE_Raw)

## adjust this value for cpm cut-off of all genes that passes filter
cpm_cutoff <- 0

## adjust this value to require the gene to appear in at least 'cpm_sample_count'
## number of samples within EITHER of the group (etoh or vehicle) to pass 
## the filter (range in [0, 6])
cpm_sample_count <- 2 # range [0, 6]

## filter out low-expressed genes according to 'cpm_cutoff' and 'cpm_sample_count'
DGE_Raw <- 
  DGE_Raw[(rowSums(cpm_gene_count[, samples$treatment %in% "etoh"] >= cpm_cutoff) >= cpm_sample_count) |
            (rowSums(cpm_gene_count[, samples$treatment %in% "vehicle"] >= cpm_cutoff) >= cpm_sample_count)
          , ]

DGE_Raw <- calcNormFactors(DGE_Raw)
DGE_Raw <- estimateDisp(DGE_Raw, 
                        design = design_samples,
                        robust = T)

## make MD plots of all 12 samples and observe the expression profile
# par(mfrow = c(3,4))
# i <- 1
# for(i in 1:12) {
#   plotMD(cpm(DGE_Raw, log = T),
#          column = i) +
#     abline(h = 0, col = "red", lty = 2, lwd = 2)
# }
# ## You can save the plot here before wiping it out in the next line
# dev.off()

## make MDS and BCV plots to observe the clustering pattern of samples
## with the cell-type element regressed out as blocking factor.
## !! TRY DIFFERENT 'cpm_cutoff' and 'cpm_sample_count' VALUES !!
## !! AND SEE THEIR EFFECTS ON the MDS PLOT.                   !!
# plotMDS(DGE_Raw)
# plotBCV(DGE_Raw)

## Run glmQLFit and glmQLFTest for DE analysis. 
## glmQLFTest is generally (but not always) more conservative. 
## Separate lines of code for glmLRTest is provided at the end 
## of the script, which is *generally* more aggressive in 
## identifying DE genes.
DGE_QLF <- glmQLFit(DGE_Raw,
                    design = design_samples,
                    robust = T)
DGE_QLF <- glmQLFTest(DGE_QLF,
                      contrast = c(1, -1, 0)) # refer to the structure of the design matrix
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
            file = "GO_results_BP_table.txt",
            quote = F, row.names = T, col.names = T,
            sep = "\t")

# GO_results_BP <- GO_results_BP[1:20, ]


## make some plots to see data distribution
plotQLDisp(DGE_QLF)
plotMD(DGE_QLF)

## export data as a data.frame
DGE_QLF_results <- DGE_QLF$table
DGE_QLF_results$Geneid <- rownames(DGE_QLF_results)
DGE_QLF_results$FDR <- p.adjust(DGE_QLF_results$PValue,
                                method = "fdr")
DGE_QLF_results <- merge(DGE_QLF_results,
                         mart_export_gencode_v33,
                         by.x = "Geneid",
                         by.y = colnames(mart_export_gencode_v33)[1])
write.table(DGE_QLF_results,
            file = "cpm_samples_1/QLFTest_results_table.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")
