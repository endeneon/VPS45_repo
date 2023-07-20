# Siwei 29 Sept 2020
# organise RNASeq data of ASD

# init
library(readr)
library(edgeR)

##
RNAseq_ASD_datExpr <- 
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/raw_data/RNAseq_ASD/RNAseq_ASD_datExpr.csv")

gene_name_list <- RNAseq_ASD_datExpr$X1
RNAseq_ASD_datExpr$X1 <- NULL
RNAseq_ASD_datExpr <- RNAseq_ASD_datExpr[, order(colnames(RNAseq_ASD_datExpr))]
rownames(RNAseq_ASD_datExpr) <- gene_name_list

RNAseq_ASD_datMeta <- 
  RNAseq_ASD_datMeta[order(RNAseq_ASD_datMeta$Dissected_Sample_ID), ]

ASD_group <- as.factor(RNAseq_ASD_datMeta$Diagnosis_)
DGE_ASD_raw <- DGEList(as.matrix(RNAseq_ASD_datExpr),
                       group = ASD_group,
                       genes = rownames(RNAseq_ASD_datExpr),
                       samples = colnames(RNAseq_ASD_datExpr),
                       remove.zeros = T)
DGE_ASD_raw <- calcNormFactors(DGE_ASD_raw)
DGE_ASD_raw <- estimateDisp(DGE_ASD_raw)

ASD_design <- model.matrix(~ 0 + ASD_group)
DGE_ASD_QLF <- glmQLFit(DGE_ASD_raw, 
                        design = ASD_design)
DGE_ASD_QLF <- glmQLFTest(DGE_ASD_QLF,
                          contrast = c(1, -1))

DGE_ASD_QLFTable <- DGE_ASD_QLF$table
DGE_ASD_QLFTable$Geneid <- rownames(DGE_ASD_QLFTable)
DGE_ASD_QLFTable$FDR <- p.adjust(DGE_ASD_QLFTable$PValue,
                                 method = "fdr")
save.image(file = "ASD_raw.RData")

save(list = "DGE_ASD_QLFTable",
     file = "ASD.RData",
     compress = T)
