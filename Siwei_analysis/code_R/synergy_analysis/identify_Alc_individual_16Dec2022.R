# Siwei 16 Dec 2022
# Identify Alc and control individuals using marker genes
# MOBP, OPALIN, UGT8, ERMN

# init
library(readr)
library(edgeR)

GSE189139_raw <-
  read_delim("TF_libraries/GSE189139_featureCounts_count_matrix.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

diff_exp_alc <-
  read_table("TF_libraries/diff.exp.alc.consumption.bmi_audit_ag_ph_cov.gene.name.txt",
             col_names = FALSE)

gencode_v35_ENSG_Genename_final <-
  read_delim("gencode.v35.ENSG.Genename.final.list",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[, c(1, 5)]
colnames(gencode_v35_ENSG_Genename_final) <-
  c("Geneid", "Gene_symbol")

gencode_v35_ENSG_Genename_final$Geneid <-
  str_split(gencode_v35_ENSG_Genename_final$Geneid,
               pattern = "\\.",
               simplify = T)[, 1]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]

df_raw <-
  merge(x = GSE189139_raw,
        y = gencode_v35_ENSG_Genename_final,
        by.x = "geneid",
        by.y = "Geneid")

df_raw <-
  df_raw[!duplicated(df_raw$Gene_symbol), ]
rownames(df_raw) <- df_raw$Gene_symbol

df_4_DGE <- df_raw
df_4_DGE$geneid <- NULL
df_4_DGE$Gene_symbol <- NULL

rownames(df_4_DGE) <- df_raw$Gene_symbol

###
df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = colnames(df_4_DGE),
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_count <-
  as.data.frame(cpm(df_4_DGE))

hist(log1p(unlist(cpm_count[rownames(cpm_count) %in% "MOBP", ])))
hist(log1p(unlist()))
hist(log1p(unlist(cpm_count[rownames(cpm_count) %in% "OPALIN", ])))
hist(log1p(unlist(cpm_count[rownames(cpm_count) %in% "UGT8", ])))
hist(log1p(unlist(cpm_count[rownames(cpm_count) %in% "ERMN", ])))

hist(log1p(unlist(cpm_count[rownames(cpm_count) %in% "C5orf17", ])),
     breaks = 100)
hist(log1p(unlist(cpm_count[rownames(cpm_count) %in% "CDHR1", ])),
     breaks = 100)
cpm_count[rownames(cpm_count) %in% "OPALIN", ]
cpm_count[rownames(cpm_count) %in% "UGT8", ]
cpm_count[rownames(cpm_count) %in% "ERMN", ]

cpm_count[rownames(cpm_count) %in%
            c("MOBP", "OPALIN", "UGT8", "ERMN"), ]
cpm_count[rownames(cpm_count) %in%
            c("C5orf17", "CDHR1"), ]


diff_exp_alc <-
  read_table("TF_libraries/diff.exp.alc.consumption.bmi_audit_ag_ph_cov.gene.name.txt",
             col_names = FALSE)

diff_exp_alc <-
  diff_exp_alc[!is.na(diff_exp_alc$X7), ]
head(diff_exp_alc)
diff_exp_alc <-
  diff_exp_alc[order(diff_exp_alc$X10), ]
diff_exp_alc[diff_exp_alc$X5 %in% c("MOBP", "OPALIN", "UGT8", "ERMN"), ]
View(head(diff_exp_alc))
