# 13 Feb 2023
# calc CPM of the 20 NGN2 RNA-seq samples

# init
library(readr)
library(readxl)

library(stringr)
library(edgeR)
library(variancePartition)

library(ggplot2)

# load data
df_raw <-
  read_delim("../STAR_output1/ReadsPerGene/output/NGN_20_RNASeq_ReadsPerGene_STAR.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

df_raw$Geneid <-
  str_split(string = df_raw$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_raw <-
  df_raw[!duplicated(df_raw$Geneid), ]

gencode_gene_name_id <-
  read_delim("/nvmefs/VPS45/R_VPS45/gencode.v35.ENSG.Genename.final.list",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
colnames(gencode_gene_name_id) <-
  c("Geneid", "CHR", "START", "END", "Gene_symbol")
gencode_gene_name_id$Geneid <-
  str_split(string = gencode_gene_name_id$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
gencode_gene_name_id <-
  gencode_gene_name_id[!duplicated(gencode_gene_name_id$Geneid), ]


df_DGE <-
  DGEList(counts = df_raw[2:ncol(df_raw)],
          samples = colnames(df_raw)[2:ncol(df_raw)],
          genes = df_raw$Geneid,
          remove.zeros = T)

cpm_DGE <-
  as.data.frame(cpm(df_DGE))

colnames(cpm_DGE) <-
  str_replace(string = colnames(cpm_DGE),
              pattern = "R21_",
              replacement = "NGN2_CD_")
cpm_DGE$Geneid <-
  df_DGE$genes$genes

cpm_DGE <-
  merge(x = cpm_DGE,
        y = gencode_gene_name_id,
        by = "Geneid")

write.table(cpm_DGE,
            file = "NGN2_20_samples_cpm.tsv",
            row.names = F, col.names = T,
            quote = F, sep = "\t")
