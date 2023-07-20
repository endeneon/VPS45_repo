# Siwei 21 Sept 2021
# Match old Table S8 with new P values

# init
library(readxl)
library(readr)

# load data
Table_S8_data <- read_excel("Table_S8_data.xlsx")
Table_S8_data <- Table_S8_data[!duplicated(Table_S8_data$Geneid), ]

VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021 <-
  read_delim("VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

merged_df <- merge(x = Table_S8_data,
                   y = VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021,
                   by = "Geneid",
                   all.x = T)
write.table(merged_df,
            file = "merged_df_for_Table_S8.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
