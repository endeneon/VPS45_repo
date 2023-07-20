# 28 Sept 2020 Siwei
# organise Geschwind et al Alz data (PFC)


# Sys.setenv(CURL_CA_BUNDLE = "/home/zhangs3/Data/Anaconda3-envs/rstudio/ssl/cacert.pem")
# source("https://bioconductor.org/biocLite.R")

library(readr)
library(stringr)
library(edgeR)

GSE44770_Alz_Geschwind_PFC <- read_delim("/nvmefs/VPS45/GSE44770_Alz_Geschwind_PFC.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
mart_EntrezID <- read_delim("mart_EntrezID.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
# mart_EntrezID <- mart_EntrezID[!is.na(mart_EntrezID$`EntrezGene transcript name ID`), ]

# use the $gene column, split by "," and make a list for indexing
gene_name_index <- str_split(GSE44770_Alz_Geschwind_PFC$gene,
                             pattern = ',', 
                             simplify = F)
colnames(ENSG_coord_gene_gencodev28) <- c("Geneid", "CHR", "START", "END", "Gene_Symbol")

GSE44770_head_col <- GSE44770_Alz_Geschwind_PFC[, 1:3]
GSE44770_head_col$Gene_Symbol <- "Nomatch"

# write the loop
i <- 1
for(i in 1:nrow(GSE44770_head_col)) {
  # print(paste("i = ", i))
  GSE_genes_to_match <- gene_name_index[[i]]
  for(j in 1:nrow(ENSG_coord_gene_gencodev28)) {
    output_gene_symbol <- match(table = GSE_genes_to_match,
                                x = ENSG_coord_gene_gencodev28$Gene_Symbol[j])
    if(sum(!is.na(output_gene_symbol)) > 0) {
      print(paste("i = ", i))
      GSE44770_head_col$Gene_Symbol[i] <- ENSG_coord_gene_gencodev28$Gene_Symbol[j]
      break
    }
  }
}
save.image()
# need to add a column of index #
GSE44770_head_col$index_no <- 1:nrow(GSE44770_head_col)
# remove all "nomatch" rows
GSE44770_head_col <- GSE44770_head_col[!(GSE44770_head_col$Gene_Symbol %in% "Nomatch"), ]
# will order by the Gene_symbol column and remove duplicates
GSE44770_head_col <- GSE44770_head_col[order(GSE44770_head_col$Gene_Symbol), ]
GSE44770_head_col <- GSE44770_head_col[!duplicated(GSE44770_head_col$Gene_Symbol), ]

# assign ENSG numbers
GSE44770_head_col <- merge(GSE44770_head_col, ENSG_coord_gene_gencodev28,
                           by = "Gene_Symbol")
GSE44770_head_col <- GSE44770_head_col[order(GSE44770_head_col$index_no), ]
GSE44770_head_col <- GSE44770_head_col[!duplicated(GSE44770_head_col$Geneid), ]


GSE44770_Alz_Geschwind_PFC_ENSG <- GSE44770_Alz_Geschwind_PFC[GSE44770_head_col$index_no, ]
rownames(GSE44770_Alz_Geschwind_PFC_ENSG) <- GSE44770_head_col$Geneid

GSE44770_Alz_Geschwind_PFC_ENSG <- GSE44770_Alz_Geschwind_PFC_ENSG[, -(1:3)]

###
GSE44770_Alz_edgeR_raw <- GSE44770_Alz_Geschwind_PFC_ENSG[, 1:230]
rownames(GSE44770_Alz_edgeR_raw) <- GSE44770_head_col$Geneid
# sum(is.na(GSE44770_Alz_edgeR_raw))
GSE44770_Alz_edgeR_raw[is.na(GSE44770_Alz_edgeR_raw)] <- 0

DGE_GSE44770_Alz <- DGEList(counts = as.matrix(GSE44770_Alz_edgeR_raw), 
                            group = as.factor(GSE44678_details$characteristics_ch2),
                            genes = rownames(GSE44770_Alz_edgeR_raw), 
                            samples = colnames(GSE44770_Alz_edgeR_raw),
                            remove.zeros = T)
group <- as.factor(GSE44678_details$characteristics_ch2)

DGE_GSE44770_Alz <- calcNormFactors(DGE_GSE44770_Alz)
DGE_GSE44770_Alz <- estimateDisp(DGE_GSE44770_Alz)

design <- model.matrix(~ 0 + group)
DGE_GSE44770_Alz_QLFTest <- glmQLFit(DGE_GSE44770_Alz, design = design)
DGE_GSE44770_Alz_QLFTest <- glmQLFTest(DGE_GSE44770_Alz_QLFTest, contrast = c(1, -1))

DGE_GSE44770_Alz_QLFTest_results <- DGE_GSE44770_Alz_QLFTest$table
DGE_GSE44770_Alz_QLFTest_results$FDR <- p.adjust(DGE_GSE44770_Alz_QLFTest_results$PValue,
                                                 method = "fdr")
DGE_GSE44770_Alz_QLFTest_results <- DGE_GSE44770_Alz_QLFTest_results[order(DGE_GSE44770_Alz_QLFTest_results$PValue), ]
DGE_GSE44770_Alz_QLFTest_results$Geneid <- rownames(DGE_GSE44770_Alz_QLFTest_results)

save.image()
saveRDS(object = DGE_GSE44770_Alz_QLFTest_results,
        file = "DGE_GSE44770_Alz_QLFTest_results.RDS")
###
sum(is.na(str_match(string = gene_name_index[[5]],
                    pattern = ENSG_coord_gene_gencodev28$Gene_Symbol[1])))

# Use ROSMap data
# Integrating Gene and Protein Reveals Perturbed
# Functional Networks in Alzheimer’s Disease

