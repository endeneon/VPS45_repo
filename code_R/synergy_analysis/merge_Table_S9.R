# make new Table S9
# Siwei 27 Jan 2023

# init
library(readxl)

load("~/NVME/VPS45/R_synergy_analysis/df_results_GG_3xKD_AA_no_sg2.RData")

Supplementary_Tables_VPS45_rev2023 <-
  read_excel("Supplementary_Tables_VPS45_rev2023.xlsx",
             sheet = "S9.KD_DEG_all", skip = 2)

gene_index <-
  Supplementary_Tables_VPS45_rev2023[, 1]

results_QLM_3xKD_vs_AA$Geneid <-
  rownames(results_QLM_3xKD_vs_AA)

KD_3_aligned <-
  merge(gene_index,
        results_QLM_3xKD_vs_AA,
        by= "Geneid",
        all.x = T)

KD_3_aligned <-
  KD_3_aligned[match(Supplementary_Tables_VPS45_rev2023$Geneid,
                     KD_3_aligned$Geneid), ]

write.table(KD_3_aligned,
            file = "KD_3_aligned_to_Table_S9.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)


KD_GG_AA_triple_KD <-
  as.data.frame(cbind(results_QLM_GG_vs_AA,
                      results_QLM_3xKD_vs_AA))
KD_GG_AA_triple_KD$Geneid <-
  rownames(KD_GG_AA_triple_KD)
KD_GG_AA_triple_KD <-
  KD_GG_AA_triple_KD[KD_GG_AA_triple_KD$FDR < 0.2, ]

write.table(KD_GG_AA_triple_KD,
            file = "CRISPR_AA_GG_vs_KDx3.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)
