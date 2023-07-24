# Siwei 07 Oct 2020
# make Upset plot

# init
library(UpSetR)

##
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"),
                   header = T, sep = ";")

# upset plot with FDR cutoff
df_upset <- data.frame(Geneid = VPS45_all_lines$Geneid, 
                       stringsAsFactors = F)

# set VPS45 as the first column
df_upset$VPS45 <- 0
df_upset$VPS45[df_upset$Geneid %in% 
                 df_upset$Geneid[VPS45_all_lines$FDR < 0.05]
               [(df_upset$Geneid[VPS45_all_lines$FDR < 0.05] 
                 %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id)]] <- 1

df_upset$SCZ <- 0
# df_upset$SCZ[df_upset$Geneid %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id] <- 1
df_upset$SCZ[df_upset$Geneid %in% 
                 DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id[DER_13_Disorder_DEX_Genes_details_DGE$SCZ.fdr < 0.05]
               [(DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id[DER_13_Disorder_DEX_Genes_details_DGE$SCZ.fdr < 0.05] 
                 %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id)]] <- 1

df_upset$ASD <- 0
# df_upset$ASD[df_upset$Geneid %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id] <- 1
df_upset$ASD[df_upset$Geneid %in% 
               DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id[DER_13_Disorder_DEX_Genes_details_DGE$ASD.fdr < 0.05]
             [(DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id[DER_13_Disorder_DEX_Genes_details_DGE$ASD.fdr < 0.05] 
               %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id)]] <- 1

df_upset$Bipolar <- 0
# df_upset$Bipolar[df_upset$Geneid %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id] <- 1
df_upset$Bipolar[df_upset$Geneid %in% 
               DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id[DER_13_Disorder_DEX_Genes_details_DGE$BD.fdr < 0.05]
             [(DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id[DER_13_Disorder_DEX_Genes_details_DGE$BD.fdr < 0.05] 
               %in% DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id)]] <- 1

df_upset$Alcohol <- 0
# df_upset$Alcohol[df_upset$Geneid %in% 
#                    DGE_Alc_QLF_results$Geneid] <- 1
df_upset$Alcohol[df_upset$Geneid %in% 
                   DGE_Alc_QLF_results$Geneid[DGE_Alc_QLF_results$FDR < 0.05]
                 [(DGE_Alc_QLF_results$Geneid[DGE_Alc_QLF_results$FDR < 0.05] 
                   %in% DGE_Alc_QLF_results$Geneid)]] <- 1

df_upset$`Crohn's Disease` <- 0
# df_upset$`Crohn's Disease`[df_upset$Geneid %in% DGE_Crohn_table$Geneid] <- 1
# [df_upset$Geneid %in% DGE_Crohn_table$Geneid]
df_upset$`Crohn's Disease`[df_upset$Geneid %in% 
                             DGE_Crohn_table$Geneid[DGE_Crohn_table$FDR < 0.05]
                           [(DGE_Crohn_table$Geneid[DGE_Crohn_table$FDR < 0.05] 
                             %in% DGE_Crohn_table$Geneid)]] <- 1


# df_upset$MDD <- 0
# # df_upset$MDD[df_upset$Geneid %in% DGE_MDD2_glmQLF_results$Geneid] <- 1
# df_upset$MDD[df_upset$Geneid %in% 
#                DGE_MDD2_glmQLF_results$Geneid[DGE_MDD2_glmQLF_results$FDR < 0.05]
#              [(DGE_MDD2_glmQLF_results$Geneid[DGE_MDD2_glmQLF_results$FDR < 0.05] 
#                %in% DGE_MDD2_glmQLF_results$Geneid)]] <- 1

colnames(df_upset)[7] <- "Crohn_s_Disease"
upset(df_upset,
      sets = c("VPS45", "SCZ", "ASD", "Bip",
               # "Crohn_s_Disease", 
               "Alcohol", "MDD"),
      # intersections = list(list("VPS45", "SCZ")),
      #                      list("VPS45", "ASD"),
      #                      list("VPS45", "Bip"),
      #                      list("VPS45", colnames(df_upset)[7]),
      #                      list("VPS45", "Alcohol"),
      #                      list("VPS45", "MDD")),
      keep.order = T,
      order.by = "freq",
      number.angles = 45,
      nintersects = 100)

save.image(file = "upsetR_plot.RData")


df_SCZ_ASD_Bip <- 
  data.frame(Geneid = DER_13_Disorder_DEX_Genes_details_DGE$ensembl_gene_id,
             SCZ = 0,
             ASD = 0,
             Bip = 0,
             stringsAsFactors = F)
df_SCZ_ASD_Bip$SCZ[DER_13_Disorder_DEX_Genes_details_DGE$SCZ.fdr < 0.05] <- 1
df_SCZ_ASD_Bip$ASD[DER_13_Disorder_DEX_Genes_details_DGE$ASD.fdr < 0.05] <- 1
df_SCZ_ASD_Bip$Bip[DER_13_Disorder_DEX_Genes_details_DGE$BD.fdr < 0.05] <- 1

df_upset <- data.frame(Geneid = VPS45_all_lines$Geneid, 
                       stringsAsFactors = F)
df_upset$VPS45 <- 0
df_upset$VPS45[VPS45_all_lines$FDR < 0.05] <- 1

df_upset <- merge(df_upset,
                  df_SCZ_ASD_Bip,
                  all.x = T,
                  # all.y = T,
                  by = "Geneid")
df_Alc <- data.frame(Geneid = DGE_Alc_QLF_results$Geneid,
                     Alcohol = 0)
df_Alc$Alcohol[DGE_Alc_QLF_results$FDR < 0.05] <- 1
df_upset <- merge(df_upset,
                  df_Alc,
                  all.x = T,
                  # all.y = T,
                  by = "Geneid")

df_Crohn <- data.frame(Geneid = DGE_Crohn_table$Geneid,
                       Crohn_s_Disease = 0)
df_Crohn$Crohn_s_Disease[DGE_Crohn_table$FDR < 0.05] <- 1
df_upset <- merge(df_upset,
                  df_Crohn,
                  all.x = T,
                  # all.y = T,
                  by = "Geneid")

# df_MDD2 <- data.frame(Geneid = DGE_MDD2_glmQLF_results$Geneid,
#                       MDD = 0)
# df_MDD2$MDD[DGE_MDD2_glmQLF_results$FDR < 0.05] <- 1
# df_upset <- merge(df_upset,
#                   df_MDD2,
#                   all.x = T,
#                   # all.y = T,
#                   by = "Geneid")


DGE_MDD2_glmLRT_results$Geneid <- rownames(DGE_MDD2_glmLRT_results)
df_MDD2 <- data.frame(Geneid = DGE_MDD2_glmQLF_results$Geneid,
                      MDD = 0)
df_MDD2$MDD[DGE_MDD2_glmLRT_results$FDR < 0.05] <- 1
df_upset <- merge(df_upset,
                  df_MDD2,
                  all.x = T,
                  # all.y = T,
                  by = "Geneid")
df_upset[is.na(df_upset)] <- 0

upset(df_upset,
      sets = c("VPS45", "SCZ", "ASD", "Bip",
               # "Crohn_s_Disease",
               "Alcohol", "MDD"),
      main.bar.color = "darkblue",
      sets.bar.color = "black",
      shade.color = "yellow",
      # matrix.color = "red",
      queries = list(list(query = intersects,
                          params = list(list("VPS45", "SCZ"),
                                        list("VPS45", "SCZ", "ASD"),
                                        list("VPS45", "SCZ", "ASD", "Bip"),
                                        list("VPS45", "SCZ", "ASD", "Bip", "MDD"),
                                        list("VPS45", "SCZ", "ASD", "Bip", "MDD", "Alcohol")),
                          colour = "red",
                          active = TRUE)),
      
      # intersections = list("VPS45", "SCZ", "ASD",
      #                      "Bip", "MDD", "Alcohol",
      #                      list("VPS45", "SCZ"),
      #                      list("VPS45", "ASD"),
      #                      list("VPS45", "Bip"),
      #                      # list("VPS45", colnames(df_upset)[7]),
      #                      list("VPS45", "Alcohol"),
      #                      list("VPS45", "MDD")),
      keep.order = F,
      decreasing = c(T, F),
      # order.by = "freq",
      number.angles = 90,
      nintersects = 100)

##
df_upset[is.na(df_upset)] <- 0
colnames(df_upset)[3] <- "SZ"
colnames(df_upset)[6] <- "AAD"
colnames(df_upset)[5] <- "BD"

upset(df_upset,
      sets = c("VPS45", "SZ", "ASD", "BD",
               # "Crohn_s_Disease",
               "AAD", "MDD"),
      main.bar.color = "darkblue",
      sets.bar.color = "black",
      shade.color = "yellow",
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.2, 1),
      line.size = 1.2,
      # matrix.color = "red",
      intersections = list("VPS45", "SZ", "ASD",
                           "BD", "MDD", "AAD",
                           list("ASD", "SZ"),
                           list("SZ", "BD"),
                           list("VPS45", "SZ"),
                           list("VPS45", "ASD"),
                           list("SZ", "AAD"),
                           list("SZ", "MDD"),
                           list("VPS45", "BD"),
                           list("VPS45", "AAD"),
                           list("ASD", "MDD"),
                           list("ASD", "BD"),
                           list("VPS45", "MDD"),
                           list("ASD", "AAD"),
                           list("MDD", "BD"),
                           list("VPS45", "SZ", "ASD"),
                           list("VPS45", "SZ", "BD"),
                           list("VPS45", "SZ", "AAD"),
                           list("VPS45", "SZ", "MDD"),
                           list("VPS45", "SZ", "ASD", "BD"),
                           list("VPS45", "SZ", "ASD", "AAD"),
                           list("VPS45", "SZ", "BD", "MDD"),
                           list("VPS45", "SZ", "ASD", "MDD")),
      queries = list(list(query = elements,
                          params = list("VPS45", "SZ"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "ASD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "BD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "AAD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "MDD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "ASD", "BD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "ASD", "AAD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "BD", "MDD"),
                          color = "red",
                          active = TRUE),
                     list(query = intersects,
                          params = list("VPS45", "SZ", "ASD", "MDD"),
                          color = "red",
                          active = TRUE)
                     ),
      decreasing = c(T, F),
      sets.x.label = "DE genes (FDR<0.05)",
      mainbar.y.label = "Overlapped gene counts",
      # order.by = "freq",
      number.angles = 0)

save.image(file = "upsetR_plot.RData")

sum(df_upset$SZ)
