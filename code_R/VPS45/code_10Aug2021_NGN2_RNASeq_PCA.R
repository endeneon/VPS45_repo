# Siwei 10 Aug 2021
# make PCA plot for NGN2 neuron based on RNA-Seq reads

# init
library(readr)

library(ggplot2)
library(gplots)

library(factoextra)
library(stringr)

library(edgeR)

# load data
df_raw_input <- 
  read_delim("NGN_20_RNASeq_ReadsPerGene_STAR.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

df_read_count <- df_raw_input
rownames(df_read_count) <- df_read_count$Geneid
df_read_count$Geneid <- NULL

df_read_count <-
  df_read_count[rowSums(df_read_count) != 0, ]

# calculate PCA use prcomp
p <- prcomp(as.matrix(t(df_read_count)),
            center = T,
            scale. = T,
            tol = 0)
scr_p <- fviz_eig(p)
fviz_pca_ind(p,
             palette = "Dark2",
             # label = "none",
             # addEllipses = T,
             invisible = "quali")


fviz_pca_var(p,
             geom.var = "text",
             label = "none")

#### use edgeR #####

df_read_count <- df_raw_input
Gene_id <- df_read_count$Geneid
rownames(df_read_count) <- Gene_id
df_read_count$Geneid <- NULL

DGE_df <- DGEList(counts = as.matrix(df_read_count), 
                  samples = colnames(df_read_count),
                  genes = Gene_id, 
                  remove.zeros = T)

DGE_df <- 
  DGE_df[rowSums(cpm(DGE_df) > 1) > 15, , keep.lib.sizes = F]

DGE_df <- calcNormFactors(DGE_df)
DGE_df <- estimateDisp(DGE_df)

DGE_df_cpm <- cpm(DGE_df)

# calculate PCA use prcomp
p <- prcomp(as.matrix(t(DGE_df_cpm)),
            center = T,
            scale. = T,
            tol = 0)
scr_p <- fviz_eig(p)
fviz_pca_ind(p,
             palette = "Dark2",
             # label = "none",
             # addEllipses = T,
             invisible = "quali")

DGE_df_cpm <- as.data.frame(DGE_df_cpm)
rownames(DGE_df_cpm) <- DGE_df$genes$genes

View(DGE_df_cpm[str_detect(rownames(DGE_df_cpm), # MAP2
                           "ENSG00000078018", 
                           negate = F), ])

View(DGE_df_cpm[str_detect(rownames(DGE_df_cpm), # MAP2
                           "ENSG00000111640", 
                           negate = F), ])

View(DGE_df_cpm[str_detect(rownames(DGE_df_cpm), 
                           "ENSG00000078018", # MAP2
                           negate = F) | 
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000104888", # vGlut1
                             negate = F) |
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000128683", # GAD1
                             negate = F) |
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000180176", # TH
                             negate = F) | 
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000091664", # vGlut2
                             negate = F) |
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000070748", # CHAT
                             negate = F) | 
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000131095", # GFAP
                             negate = F) | 
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000139287", # TPH2
                             negate = F) | 
                  str_detect(rownames(DGE_df_cpm), 
                             "ENSG00000205927", # OLIG2
                             negate = F) , ])

df_gene_panel <-
  DGE_df_cpm[str_detect(rownames(DGE_df_cpm), 
                        "ENSG00000078018", # MAP2
                        negate = F) | 
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000104888", # vGlut1
                          negate = F) |
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000128683", # GAD1
                          negate = F) |
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000180176", # TH
                          negate = F) | 
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000091664", # vGlut2
                          negate = F) |
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000070748", # CHAT
                          negate = F) | 
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000131095", # GFAP
                          negate = F) | 
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000139287", # TPH2
                          negate = F) | 
               str_detect(rownames(DGE_df_cpm), 
                          "ENSG00000205927", # OLIG2
                          negate = F) , ]

write.table(df_gene_panel,
            file = "NGN2_Glut_gene_panel.txt",
            quote = F, sep = "\t",
            col.names = T, row.names = T,)
