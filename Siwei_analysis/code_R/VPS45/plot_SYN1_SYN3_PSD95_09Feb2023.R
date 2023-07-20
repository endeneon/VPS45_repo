# Siwei 09 Feb 2023
# plot the expression of SYN1, SYN3, and PSD-95 in AA/GG bulk RNA-seq

# init
library(edgeR)
library(readr)

library(Rfast)
library(factoextra)
library(stringr)
library(MASS)

library(magrittr)
library(dplyr)
library(ggpubr)

library(ggplot2)
library(gplots)
library(RColorBrewer)

library(sva)
library(FactoMineR)
library(future)

library(corrplot)
library(ggrepel)
library(reshape2)
library(pheatmap)


# load AA/GG data from single KD RNA-seq

df_raw_VPS45_single_KD <- 
  read_delim("ReadsPerGene_STAR_full.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)

# add gene symbols
gencode_v35_ENSG_Genename_final <- 
  read_delim("gencode.v35.ENSG.Genename.final.list", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[, c(1, 5)]
colnames(gencode_v35_ENSG_Genename_final) <-
  c("Geneid", "Gene_symbol")

gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]

df_master <-
  df_raw_VPS45_single_KD[, c(1, 5:7, 18:20, 27:29)]

df_master <-
  merge(x = df_master,
        y = gencode_v35_ENSG_Genename_final,
        by = "Geneid")
df_master <-
  df_master[!duplicated(df_master$Gene_symbol), ]

geneid_list <- df_master$Geneid
genesymbol_list <- df_master$Gene_symbol

df_master$Geneid <- NULL
df_master$Gene_symbol <- NULL

rownames(df_master) <- genesymbol_list

df_metadata <-
  data.frame(genotype = c(rep_len("AA", length.out = 6),
                          rep_len("GG", length.out = 3)))

DGE_raw <-
  DGEList(counts = as.matrix(df_master),
          samples = colnames(df_master),
          group = df_metadata$genotype,
          genes = rownames(df_master),
          remove.zeros = T)

cpm_gene_count <- cpm(DGE_raw)

cpm_to_plot <-
  as.data.frame(cpm_gene_count[rownames(cpm_gene_count) %in% c("SYN1", 
                                                               "SYN3",
                                                               "DLG4"), ])

# 1=SYN1, 2=DLG4, 3=SYN3
df_to_plot <-
  data.frame(CPM = unlist(cpm_to_plot[1, ]),
             genotype = factor(c("AA", "AA", "AA",
                                 "AA", "AA", "AA",
                                 "GG", "GG", "GG")))
t.test(x = df_to_plot$CPM[1:6],
       y = df_to_plot$CPM[7:9], 
       alternative = "t", 
       var.equal = F)

df_to_plot <-
  data.frame(CPM = unlist(cpm_to_plot[2, ]),
             genotype = factor(c("AA", "AA", "AA",
                                 "AA", "AA", "AA",
                                 "GG", "GG", "GG")))
t.test(x = df_to_plot$CPM[1:6],
       y = df_to_plot$CPM[7:9], 
       alternative = "t", 
       var.equal = F)

df_to_plot <-
  data.frame(CPM = unlist(cpm_to_plot[3, ]),
             genotype = factor(c("AA", "AA", "AA",
                                 "AA", "AA", "AA",
                                 "GG", "GG", "GG")))
t.test(x = df_to_plot$CPM[1:6],
       y = df_to_plot$CPM[7:9], 
       alternative = "t", 
       var.equal = F)

ggplot(df_to_plot,
       aes(x = genotype,
           y = CPM,
           fill = genotype,
           group = genotype)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(size = 1,
              width = -1) +
  ylim(0, max(df_to_plot$CPM) + 10) +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")[c(1, 3)]) +
  theme_classic() +
  theme(axis.text = element_text(size = 12))
