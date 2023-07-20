# 29 Sept 2020
# organise merged GTEx+PsychEncode RNA-Seq data

# init
library(stringr)
library(edgeR)


gene_name_index <- DER_01_PEC_GTEx_merged_Gene_expression_matrix_normalized$gene_id
gene_name_index <- str_split(gene_name_index, 
                             pattern = "\\.",
                             simplify = T)[, 1]
# gene_name_index <- as.data.frame(gene_name_index)
DER_01_PEC_GTEx_merged_Gene_expression_matrix_normalized$gene_id <- NULL


PEG_GTEx_phenotypes_metadata_matched <- 
  PEG_GTEx_phenotypes_metadata[match(colnames(DER_01_PEC_GTEx_merged_Gene_expression_matrix_normalized),
                                     PEG_GTEx_phenotypes_metadata$individualID), ]

Exp_matrix <- 
  DER_01_PEC_GTEx_merged_Gene_expression_matrix_normalized[, !(is.na(PEG_GTEx_phenotypes_metadata_matched$individualID))]
Meta_matrix <- PEG_GTEx_phenotypes_metadata_matched[!(is.na(PEG_GTEx_phenotypes_metadata_matched$individualID)), ]

##
Exp_matrix <- Exp_matrix[, order(colnames(Exp_matrix))]
Meta_matrix <- Meta_matrix[order(Meta_matrix$individualID), ]

Meta_matrix_CMC_excluded <- 
  Meta_matrix[!(Meta_matrix$individualID %in% CMC_MSSM_Penn_Pitt_Clinical$Individual_ID), ]
Exp_matrix_CMC_excluded <- 
  Exp_matrix[, colnames(Exp_matrix) %in% Meta_matrix_CMC_excluded$individualID]

Meta_matrix_CMC_excluded <- 
  Meta_matrix_CMC_excluded[Meta_matrix_CMC_excluded$diagnosis %in% c("Control",
                                                                     "Schizophrenia")
                           , ]
Exp_matrix_CMC_excluded <- 
  Exp_matrix_CMC_excluded[, colnames(Exp_matrix_CMC_excluded) %in% Meta_matrix_CMC_excluded$individualID]

Meta_matrix_CMC_excluded <- 
  Meta_matrix_CMC_excluded[order(Meta_matrix_CMC_excluded$individualID), ]
Exp_matrix_CMC_excluded <-
  Exp_matrix_CMC_excluded[, order(colnames(Exp_matrix_CMC_excluded))]
rownames(Exp_matrix_CMC_excluded) <- gene_name_index

group_exp <- as.factor(Meta_matrix_CMC_excluded$diagnosis)

## Need to use the native limma and eBayes() functions here 
## since the counts were normalised FPKM already
## convert FPKM back to pseudocounts
exp_design <- model.matrix(~ 0 + group_exp)

limma_matrix_CMC_excluded <- lmFit(Exp_matrix_CMC_excluded,
                                   design = exp_design)
limma_matrix_CMC_excluded <- eBayes(limma_matrix_CMC_excluded)

CMC_excluded_output <- topTable(limma_matrix_CMC_excluded,
                                coef = 2, 
                                number = nrow(Exp_matrix_CMC_excluded))
CMC_excluded_output$Geneid <- rownames(CMC_excluded_output)
saveRDS(CMC_excluded_output,
        file = "CMC_excluded_output.RDS")


##
save.image(file = "GTEx_PsychENCODE_raw.RData")
