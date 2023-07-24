# Siwei 29 Sept 2020
# test correlations

# init


## raw import data used AA as baseline, transformation needed
all_gene_list_line_11$logFC <- 0 - all_gene_list_line_11$logFC
all_gene_list_line_12$logFC <- 0 - all_gene_list_line_12$logFC



DGE_CMC_Alz_GSE44770_glmQLFtest <- readRDS("DGE_GSE44770_Alz_QLFTest_results.RDS")

## test correlation, use alternative = "two.sided"

### test use VPS45 CRISPR edited line, GG as baseline
# SCZ_CMC
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_CMC_SCZ_glmQLFTest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

test_corr_var <- merge(VPS45_all_lines, 
                       CMC_DGE_QLFTable,
                       by = "Geneid")
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)

# Alz
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_CMC_Alz_GSE44770_glmQLFtest,
                       by = "Geneid")
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG

test_corr_var <- test_corr_var[test_corr_var$FDR.x < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$FDR.y < 0.05, ]

cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# BP
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_CMC_BP_glmQLFTest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# AFF
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_CMC_AFF_glmQLFTest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# GVEx_Bipolar
test_corr_var <- merge(VPS45_all_lines, 
                       GVEx_Bipolar_glmQLFTest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# GVEx_SCZ
test_corr_var <- merge(VPS45_all_lines, 
                       GVEx_SCZ_glmQLFTest,
                       by = "Geneid")
test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# CMC_excluded_SCZ
SCZ_CMC_excluded <- readRDS(file = "CMC_excluded_output.RDS")
SCZ_CMC_excluded$logFC <- 0 - SCZ_CMC_excluded$logFC
test_corr_var <- merge(VPS45_all_lines, 
                       SCZ_CMC_excluded,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# ASD
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_ASD_QLFTable,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)


## Alz RNASeq Cell Reports
test_corr_var <- merge(VPS45_all_lines, 
                       X1_s2_0_S2211124719308538_mmc2,
                       by.x = "Geneid",
                       by.y = "Ensembl ID")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$log2.AD.NDC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)


## AAD Microarray

test_corr_var <- merge(VPS45_all_lines, 
                       Microarray_AAD_metaanalysis_092017,
                       by.x = "Geneid",
                       by.y = "X1")

test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- test_corr_var[test_corr_var$FDR < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$fdr < 0.05, ]


cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$SE,
         alternative = "t",
         method = "s",
         exact = F)


0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$beta,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## AAD RNASeq Kapoor et al 2019

test_corr_var <- merge(VPS45_all_lines, 
                       Alc_vs_ctrl,
                       by.x = "Geneid",
                       by.y = "id")

test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- test_corr_var[test_corr_var$FDR < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$padj < 0.05, ]


cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$log2FoldChange,
         alternative = "t",
         method = "s",
         exact = F)


0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$log2FoldChange,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## AAD RNASeq David Goldman 2014

test_corr_var <- merge(VPS45_all_lines, 
                       DGE_Alc_QLF_results,
                       by.x = "Geneid",
                       by.y = "Geneid")

test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- test_corr_var[test_corr_var$FDR.x < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$FDR.y < 0.05, ]


cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)


0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$log2FoldChange,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## MDD Microarray

test_corr_var <- merge(VPS45_all_lines, 
                       Microarray_MDD_metaanalysis_092017,
                       by.x = "Geneid",
                       by.y = "X1")

test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- test_corr_var[test_corr_var$FDR < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$fdr < 0.05, ]


cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$SE,
         alternative = "t",
         method = "s",
         exact = F)


0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$beta * 10,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## MDD microarray use customised analysis
test_corr_var <- merge(VPS45_all_lines, 
                       exp_MDD_MDD,
                       by.x = "Geneid",
                       by.y = "Geneid")


test_corr_var$fdr <- p.adjust(test_corr_var$P.Value, 
                              method = "fdr")


test_corr_var <- test_corr_var[test_corr_var$FDR < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$fdr < 0.05, ]


cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)


0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)



## IBD Microarray

test_corr_var <- merge(VPS45_all_lines, 
                       Microarray_IBD_metaanalysis_092017,
                       by.x = "Geneid",
                       by.y = "X1")

test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- test_corr_var[test_corr_var$FDR < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$fdr < 0.05, ]


cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$beta,
         alternative = "t",
         method = "s",
         exact = F)


0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$beta,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## Crohn's RNASeq
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_Crohn_table,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG

test_corr_var <- test_corr_var[test_corr_var$FDR.x < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$FDR.y < 0.05, ]

cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## MDD RNASeq
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_MDD_QLFResults,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG

test_corr_var <- test_corr_var[test_corr_var$FDR.x < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$FDR.y < 0.05, ]

cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## MDD RNASeq glmLRT
DGE_MDD2_glmLRT_results$Geneid <- rownames(DGE_MDD2_glmLRT_results)
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_MDD2_glmLRT_results,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG

test_corr_var <- test_corr_var[test_corr_var$FDR.x < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$PValue.y < 0.05, ]

cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "g",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

## MDD2 RNASeq
test_corr_var <- merge(VPS45_all_lines, 
                       DGE_MDD2_glmQLF_results,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG

test_corr_var <- test_corr_var[test_corr_var$FDR.x < 0.05, ]
test_corr_var <- test_corr_var[test_corr_var$FDR.y < 1, ]

cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)



######


### test use VPS45 CROPSeq KD, GG as baseline


# SCZ_CMC
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       DGE_CMC_SCZ_glmQLFTest,
                       by = "Geneid")

# test_corr_var <- test_corr_var[test_corr_var$FDR.y < 0.05, ]
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       CMC_DGE_QLFTable,
                       by = "Geneid")

# test_corr_var <- test_corr_var[test_corr_var$FDR.y < 0.05, ]
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# Alz
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       DGE_CMC_Alz_GSE44770_glmQLFtest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# BP
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       GVEx_Bipolar_glmQLFTest,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.05, ]
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# AFF
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       DGE_CMC_AFF_glmQLFTest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# GTEx_Bipolar
test_corr_var <- merge(VPS45_all_lines, 
                       GVEx_Bipolar_glmQLFTest,
                       by = "Geneid")
cor.test(test_corr_var$logFC.groupGG,
         test_corr_var$logFC,
         alternative = "t",
         method = "s",
         exact = F)
0 - log10(cor.test(test_corr_var$logFC.groupGG,
                   test_corr_var$logFC,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# GTEx_SCZ
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       GVEx_SCZ_glmQLFTest,
                       by = "Geneid")
test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# CMC_excluded_SCZ
SCZ_CMC_excluded <- readRDS(file = "CMC_excluded_output.RDS")
SCZ_CMC_excluded$logFC <- 0 - SCZ_CMC_excluded$logFC
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       SCZ_CMC_excluded,
                       by = "Geneid")
test_corr_var <- test_corr_var[test_corr_var$PValue < 0.05, ]
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

# ASD
test_corr_var <- merge(gRNA_neuron_CROPSeq_VPS45_2_gene_table, 
                       DGE_ASD_QLFTable,
                       by = "Geneid")
# test_corr_var <- test_corr_var[test_corr_var$PValue.x < 0.01, ]
cor.test(test_corr_var$logFC.x,
         test_corr_var$logFC.y,
         alternative = "t",
         method = "s",
         exact = F)

0 - log10(cor.test(test_corr_var$logFC.x,
                   test_corr_var$logFC.y,
                   alternative = "t",
                   method = "s",
                   exact = F)$p.value)

