# Siwei 08 Dec 2022
# quick test for count tables

# init
library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)

library(ggplot2)
library(gplots)
library(RColorBrewer)


gene_count_table <- 
  read_delim("ReadsPerGene_STAR_full.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)

gene_count_table <- 
  read_delim("VPS45_3xKD_ReadsPerGene_STAR_08Dec2022_stranded_correct.txt", 
             "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE)


df_raw <-
  gene_count_table

df_raw$Geneid <-
  str_split(df_raw$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]

df_raw <-
  df_raw[!duplicated(df_raw$Geneid), ]

df_4_DGE <- df_raw
# df_4_DGE$Geneid

df_4_DGE$Geneid <- NULL
rownames(df_4_DGE) <- df_raw$Geneid

### make DGEList
df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = colnames(df_4_DGE),
          genes = rownames(df_4_DGE),
          remove.zeros = T)

DGE_cpm <-
  as.data.frame(cpm(df_4_DGE))

DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000285184", ] # lncRNA
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000136631", ] # VPS45
DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000118292", ] # C1orf54

DGE_cpm[rownames(DGE_cpm) %in% "ENSG00000284202", ] 

##
View(gene_count_table[str_detect(string = gene_count_table$Geneid,
                                 pattern = "ENSG00000118292"), ])







#####
df_raw <-
  read_table("FRiP_results.txt", 
             col_names = FALSE)
df_raw$X3 <- NULL

colnames(df_raw) <-
  c("Sample", "Total_Reads", "Reads_in_Peaks")

df_raw$Sample <-
  str_split(string = df_raw$Sample,
            pattern = "_",
            simplify = T)[, 1]

df_raw$Fraction_Reads_in_Peaks <-
  df_raw$Reads_in_Peaks / df_raw$Total_Reads

df_raw$Sample <-
  str_replace(string = df_raw$Sample,
              pattern = "Ast",
              replacement = "Astrocyte")
df_raw$Sample <-
  str_replace(string = df_raw$Sample,
              pattern = "GA",
              replacement = "GABAergic_neuron")
df_raw$Sample <-
  str_replace(string = df_raw$Sample,
              pattern = "MG",
              replacement = "Microglia")
df_raw$Sample <-
  str_replace(string = df_raw$Sample,
              pattern = "R21",
              replacement = "NGN2_glut_neuron")

df_raw <-
  df_raw[order(df_raw$Sample), ]

write.table(df_raw,
            file = "FRiP_tables_bulk_ATAC.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)
