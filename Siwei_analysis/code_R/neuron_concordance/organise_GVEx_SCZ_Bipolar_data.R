# Siwei 29 Sept 2020
# organise RNASeq data from GVEx SCZ and Bipolar

# init
library(readr)
library(edgeR)

##
RNAseq_SCZ_BD_GVEX_datExpr <- 
  read_csv("/nvmefs/VPS45/shared_molecular_neuropathology_mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/raw_data/RNAseq_GVEX/RNAseq_SCZ_BD_GVEX_datExpr.csv")
RNAseq_SCZ_BD_GVEX_datExpr <- 
  RNAseq_SCZ_BD_GVEX_datExpr[!duplicated(RNAseq_SCZ_BD_GVEX_datExpr$X1), ]
gene_name_list <- RNAseq_SCZ_BD_GVEX_datExpr$X1

RNAseq_SCZ_BD_GVEX_datExpr$X1 <- NULL
rownames(RNAseq_SCZ_BD_GVEX_datExpr) <- gene_name_list

RNAseq_SCZ_BD_GVEX_datExpr <- 
  RNAseq_SCZ_BD_GVEX_datExpr[, order(colnames(RNAseq_SCZ_BD_GVEX_datExpr))]
rownames(RNAseq_SCZ_BD_GVEX_datExpr) <- gene_name_list

RNAseq_SCZ_BD_GVEX_datMeta <- 
  RNAseq_SCZ_BD_GVEX_datMeta[order(RNAseq_SCZ_BD_GVEX_datMeta$Individual_ID..RNAseq.library.BID.), ]

##
GVEx_group <- as.factor(RNAseq_SCZ_BD_GVEX_datMeta$Diagnosis)
DGE_GVEx_SCZ_Bipolar <- DGEList(as.matrix(RNAseq_SCZ_BD_GVEX_datExpr),
                                group = GVEx_group,
                                genes = rownames(RNAseq_SCZ_BD_GVEX_datExpr),
                                samples = colnames(RNAseq_SCZ_BD_GVEX_datExpr),
                                remove.zeros = T)

DGE_GVEx_SCZ_Bipolar <- calcNormFactors(DGE_GVEx_SCZ_Bipolar)
DGE_GVEx_SCZ_Bipolar <- estimateDisp(DGE_GVEx_SCZ_Bipolar)

GVEx_design <- model.matrix(~ 0 + GVEx_group)
DGE_GVEx_SCZ_Bipolar_glmQLFit <- glmQLFit(DGE_GVEx_SCZ_Bipolar, 
                                          design = GVEx_design)

# calculate two diseases separately
# > head(GVEx_design)
# GVEx_groupBP GVEx_groupControl GVEx_groupSCZ
# 1            1                 0             0
# 2            0                 0             1

## GVEX_SCZ
GVEx_SCZ_glmQLFTest <- glmQLFTest(DGE_GVEx_SCZ_Bipolar_glmQLFit,
                                  contrast = c(0, -1, 1))
GVEx_SCZ_glmQLFTest <- GVEx_SCZ_glmQLFTest$table
GVEx_SCZ_glmQLFTest$Geneid <- rownames(GVEx_SCZ_glmQLFTest)
GVEx_SCZ_glmQLFTest$FDR <- p.adjust(GVEx_SCZ_glmQLFTest$PValue)

## GVEX_Bipolar
GVEx_Bipolar_glmQLFTest <- glmQLFTest(DGE_GVEx_SCZ_Bipolar_glmQLFit,
                                      contrast = c(1, -1, 0))
GVEx_Bipolar_glmQLFTest <- GVEx_Bipolar_glmQLFTest$table
GVEx_Bipolar_glmQLFTest$Geneid <- rownames(GVEx_Bipolar_glmQLFTest)
GVEx_Bipolar_glmQLFTest$FDR <- p.adjust(GVEx_Bipolar_glmQLFTest$PValue)

save.image(file = "GVEx_SCZ_Bipolar_raw.RData")

save(list = c("GVEx_Bipolar_glmQLFTest",
              "GVEx_SCZ_glmQLFTest"),
     file = "GVEx_SCZ_Bipolar.RData",
     compress = T)
