# Siwei 29 Sept 2020
# Organise CMC_SCZ data

# init
library(readr)
library(edgeR)

##

CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw <- 
  read_delim("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/raw_data/RNAseq_CMC/CommonMind-release1/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv",
             "\t", escape_double = FALSE, trim_ws = TRUE)


rownames(CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw) <- 
  as.character(CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw$Geneid)
gene_names_list_raw <- rownames(CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw)
CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw <- 
  CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw[, -1]
rownames(CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw) <- gene_names_list_raw

CMC_MSSM_Penn_Pitt_Clinical_ordered <-
  CMC_MSSM_Penn_Pitt_Clinical[CMC_MSSM_Penn_Pitt_Clinical$DLPFC_RNA_Sequencing_Sample_ID %in%
                                colnames(CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw), ]
CMC_MSSM_Penn_Pitt_DLPFC_Expression <- 
  CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw[,
                                                                    colnames(CMC_MSSM_Penn_Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw) %in%
                                                                      CMC_MSSM_Penn_Pitt_Clinical_ordered$DLPFC_RNA_Sequencing_Sample_ID]
rownames(CMC_MSSM_Penn_Pitt_DLPFC_Expression) <- gene_names_list_raw

CMC_MSSM_Penn_Pitt_Clinical_ordered <- 
  CMC_MSSM_Penn_Pitt_Clinical_ordered[order(CMC_MSSM_Penn_Pitt_Clinical_ordered$DLPFC_RNA_Sequencing_Sample_ID), ]
CMC_MSSM_Penn_Pitt_DLPFC_Expression <- 
  CMC_MSSM_Penn_Pitt_DLPFC_Expression[, 
                                      order(colnames(CMC_MSSM_Penn_Pitt_DLPFC_Expression))]
rownames(CMC_MSSM_Penn_Pitt_DLPFC_Expression) <- gene_names_list_raw
###

CMC_group <- factor(CMC_MSSM_Penn_Pitt_Clinical_ordered$Dx)
DGE_CMC_DLPFC_Raw <- DGEList(as.matrix(CMC_MSSM_Penn_Pitt_DLPFC_Expression),
                             group = CMC_group,
                             genes = rownames(CMC_MSSM_Penn_Pitt_DLPFC_Expression),
                             samples = colnames(CMC_MSSM_Penn_Pitt_DLPFC_Expression),
                             remove.zeros = T)
DGE_CMC_DLPFC_Raw <- calcNormFactors(DGE_CMC_DLPFC_Raw)
DGE_CMC_DLPFC_Raw <- estimateDisp(DGE_CMC_DLPFC_Raw)

CMC_design <- model.matrix(~ 0 + CMC_group)
DGE_CMC_DLPFC_QLFit <- glmQLFit(DGE_CMC_DLPFC_Raw,
                                design = CMC_design)

# calculate different diseases, control position = 3
# > head(CMC_design)
# CMC_groupAFF CMC_groupBP CMC_groupControl CMC_groupSCZ
# 1            0           1                0            0
# 2            0           1                0            0
# 3            1           0                0            0

## SCZ
DGE_CMC_SCZ_glmQLFTest <- glmQLFTest(DGE_CMC_DLPFC_QLFit,
                                     contrast = c(0, 0, -1, 1))
DGE_CMC_SCZ_glmQLFTest <- DGE_CMC_SCZ_glmQLFTest$table
DGE_CMC_SCZ_glmQLFTest$Geneid <- rownames(DGE_CMC_SCZ_glmQLFTest)
DGE_CMC_SCZ_glmQLFTest$FDR <- p.adjust(DGE_CMC_SCZ_glmQLFTest$PValue,
                                       method = "fdr")

## BP
DGE_CMC_BP_glmQLFTest <- glmQLFTest(DGE_CMC_DLPFC_QLFit,
                                     contrast = c(0, 1, -1, 0))
DGE_CMC_BP_glmQLFTest <- DGE_CMC_BP_glmQLFTest$table
DGE_CMC_BP_glmQLFTest$Geneid <- rownames(DGE_CMC_BP_glmQLFTest)
DGE_CMC_BP_glmQLFTest$FDR <- p.adjust(DGE_CMC_BP_glmQLFTest$PValue,
                                       method = "fdr")

## AFF
DGE_CMC_AFF_glmQLFTest <- glmQLFTest(DGE_CMC_DLPFC_QLFit,
                                    contrast = c(1, 0, -1, 0))
DGE_CMC_AFF_glmQLFTest <- DGE_CMC_AFF_glmQLFTest$table
DGE_CMC_AFF_glmQLFTest$Geneid <- rownames(DGE_CMC_AFF_glmQLFTest)
DGE_CMC_AFF_glmQLFTest$FDR <- p.adjust(DGE_CMC_AFF_glmQLFTest$PValue,
                                      method = "fdr")


save.image(file = "CMC_SCZ.RData")

save(list = c("DGE_CMC_AFF_glmQLFTest",
              "DGE_CMC_BP_glmQLFTest",
              "DGE_CMC_SCZ_glmQLFTest"),
           file = "CMC_SCZ_AFF_Bipolar.RData",
           compress = T)
