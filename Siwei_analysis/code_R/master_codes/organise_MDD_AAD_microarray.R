# Siwei 01 Oct 2020
# organize MDD and AAD microarray data
# note microarray data has beta not logFC, 

# init
library(edgeR)

# remove NAs
Microarray_AAD_metaanalysis_092017 <- 
  Microarray_AAD_metaanalysis_092017[!(is.na(Microarray_AAD_metaanalysis_092017$beta)), ]

sum(Microarray_AAD_metaanalysis_092017$beta)


sum(is.na(Microarray_AAD_metaanalysis_092017))

# MDD data use the .RData directly from the database
# Sibillel et al.

MDD_group <- as.factor(datMeta$GROUP)
MDD_design <- model.matrix(~ 0 + MDD_group)

lmFit_MDD <- lmFit(datExpr, 
                   design = MDD_design)
lmFit_MDD <- eBayes(lmFit_MDD)
exp_MDD_ctrl <- topTable(lmFit_MDD,
                         coef = 1,
                         number = nrow(lmFit_MDD))
exp_MDD_MDD <- topTable(lmFit_MDD,
                        coef = 2,
                        number = nrow(lmFit_MDD))

exp_MDD_MDD$`-logFC` <- 0 - (exp_MDD_MDD$logFC-exp_MDD_ctrl$logFC)
exp_MDD_MDD$Geneid <- rownames(exp_MDD_MDD)

####

test_corr_var <- merge(Microarray_MDD_metaanalysis_092017, 
                       Microarray_SCZ_metaanalysis_092017,
                       by = "X1")
test_corr_var <- merge(Microarray_ASD_metaanalysis_092017, 
                       Microarray_SCZ_metaanalysis_092017,
                       by = "X1")
test_corr_var <- merge(Microarray_BD_metaanalysis_092017, 
                       Microarray_SCZ_metaanalysis_092017,
                       by = "X1")
test_corr_var <- merge(Microarray_IBD_metaanalysis_092017, 
                       Microarray_SCZ_metaanalysis_092017,
                       by = "X1")

# test_corr_var$fdr <- p.adjust(test_corr_var$P.Value, 
#                               method = "fdr")
# 
# 
# test_corr_var <- test_corr_var[test_corr_var$FDR < 0.05, ]
# test_corr_var <- test_corr_var[test_corr_var$fdr < 0.05, ]


cor.test(test_corr_var$beta.x,
         test_corr_var$beta.y,
         alternative = "t",
         method = "s",
         exact = F)


##
library(Rfast)

exp_MDD_matrix <- datExpr[, datMeta$GROUP %in% "MDD"]
exp_ctrl_matrix <- datExpr[, datMeta$GROUP %in% "Control"]

exp_cor_matrix <- data.frame(ctrl = rep(0, nrow(datExpr)))
exp_cor_matrix$ctrl <- rowMeans(exp_ctrl_matrix)
exp_cor_matrix$MDD <- rowMeans(exp_MDD_matrix)
exp_cor_matrix$`-logFC` <- exp_cor_matrix$MDD - exp_cor_matrix$ctrl

exp_cor_matrix$Geneid <- rownames(datExpr)
exp_cor_matrix <- 
  exp_cor_matrix[exp_cor_matrix$Geneid %in% 
                   Microarray_MDD_metaanalysis_092017$X1
                 , ]

test_corr_var <- merge(VPS45_all_lines,
                       exp_cor_matrix,
                       by = "Geneid")
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG

cor.test(test_corr_var$`-logFC`,
         test_corr_var$logFC.groupGG,
         alternative = "t",
         method = "s",
         exact = F)

Microarray_MDD_metaanalysis_092017 <- 
  Microarray_MDD_metaanalysis_092017[Microarray_MDD_metaanalysis_092017$fdr < 0.05, ]



## AAD
exp_EtOH_matrix <- datExpr[, datMeta$Group %in% "ETOH"]
exp_ctrl_matrix <- datExpr[, datMeta$Group %in% "CTL"]

exp_cor_matrix <- data.frame(ctrl = rep(0, nrow(datExpr)))
exp_cor_matrix$ctrl <- rowMeans(exp_ctrl_matrix)
exp_cor_matrix$EtOH <- rowMeans(exp_EtOH_matrix)
exp_cor_matrix$`-logFC` <- exp_cor_matrix$EtOH - exp_cor_matrix$ctrl

exp_cor_matrix$Geneid <- rownames(datExpr)
exp_cor_matrix <- 
  exp_cor_matrix[exp_cor_matrix$Geneid %in%
                   Microarray_AAD_metaanalysis_092017$X1
                 , ]



test_corr_var <- merge(VPS45_all_lines,
                       exp_cor_matrix,
                       by = "Geneid")
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- 
  test_corr_var[test_corr_var$FDR < 0.05, ]


cor.test(test_corr_var$`-logFC`,
         test_corr_var$logFC.groupGG,
         alternative = "t",
         method = "s",
         exact = F)

Microarray_AAD_metaanalysis_092017 <- 
  Microarray_AAD_metaanalysis_092017[Microarray_AAD_metaanalysis_092017$fdr < 0.05, ]


## IBD
exp_IBD_matrix <- datExpr[, datMeta$Characteristics.DiseaseState. %in% "diseased"]
exp_ctrl_matrix <- datExpr[, datMeta$Characteristics.DiseaseState. %in% "unaffected"]




exp_cor_matrix <- data.frame(ctrl = rep(0, nrow(datExpr)))
exp_cor_matrix$ctrl <- rowMeans(exp_ctrl_matrix)
exp_cor_matrix$IBD <- rowMeans(exp_IBD_matrix)
exp_cor_matrix$`-logFC` <- exp_cor_matrix$IBD - exp_cor_matrix$ctrl

exp_cor_matrix$Geneid <- rownames(datExpr)

test_corr_var <- merge(VPS45_all_lines,
                       exp_cor_matrix,
                       by = "Geneid")
test_corr_var$logFC.groupGG <- 0 - test_corr_var$logFC.groupGG
test_corr_var <- 
  test_corr_var[test_corr_var$FDR < 0.05, ]


cor.test(test_corr_var$`-logFC`,
         test_corr_var$logFC.groupGG,
         alternative = "t",
         method = "s",
         exact = F)

Microarray_IBD_metaanalysis_092017 <- 
  Microarray_IBD_metaanalysis_092017[Microarray_IBD_metaanalysis_092017$fdr < 0.05, ]
exp_cor_matrix <- 
  exp_cor_matrix[exp_cor_matrix$Geneid %in% 
                   Microarray_IBD_metaanalysis_092017$X1
                 ,]
