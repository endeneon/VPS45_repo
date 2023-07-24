# Siwei 14 Mar 2022
# make .gene_loc files

# init
library(stringr)
library(readr)

df_NCBI_37.3 <-
  read_delim("~/NVME/VPS45/organoids_SETD1A_hg19/MAGMA_annotation/NCBI37.3/Rev.NCBI37.3.gene.loc", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
df_NCBI_37.3$X6 <- NULL
colnames(df_NCBI_37.3) <-
  c("Gene_Symbol", "CHR", "START", "END", "STRAND")

gene_set_list <-
  list.files(path = ".",
             pattern = ".GeneSet$")

i <- 1L
for (i in 1:length(gene_set_list)) {
  current_gene_set <-
    read_delim(gene_set_list[i], 
               delim = "\t", escape_double = FALSE, 
               col_names = FALSE, trim_ws = TRUE)
  current_gene_set <- unlist(current_gene_set)
  df_filtered_gene_loc <-
    df_NCBI_37.3[df_NCBI_37.3$Gene_Symbol %in% current_gene_set, ]
  write.table(df_filtered_gene_loc,
              file = paste("gene_loc/",
                           str_replace_all(gene_set_list[i],
                                           pattern = "GeneSet",
                                           replacement = "gene_loc"),
                           sep = ""),
              quote = F, sep = "\t",
              col.names = F, row.names = F)
}

#######             
i <- 1L
for (i in 1:length(gene_set_list)) {
  current_gene_set <-
    read_delim(gene_set_list[i], 
               delim = "\t", escape_double = FALSE, 
               col_names = FALSE, trim_ws = TRUE)
  current_gene_set <- unlist(current_gene_set)
  if (i == 1) {
    all_gene_list <- current_gene_set
  } else {
    all_gene_list <- c(all_gene_list, current_gene_set)
  }
}

all_gene_list <- unique(all_gene_list)

NCBI_37.3_filtered_all_gene_list <-
  df_NCBI_37.3[df_NCBI_37.3$Gene_Symbol %in% all_gene_list, ]

write.table(NCBI_37.3_filtered_all_gene_list,
            file = "gene_loc/full_set_DE_gene_list.gene_loc",
            quote = F, sep = "\t",
            row.names = F, col.names = F)

######
i <- 1L

for (i in 1:length(gene_set_list)) {
  print(i)
  current_gene_set <-
    read_delim(gene_set_list[i], 
               delim = "\t", escape_double = FALSE, 
               col_names = FALSE, trim_ws = TRUE)
  current_gene_set <- unlist(current_gene_set)
  current_gene_set <- 
    current_gene_set[current_gene_set %in% df_NCBI_37.3$Gene_Symbol]
  
  if (i == 1) {
    set_annot_df <-
      data.frame(type = as.character(gene_set_list[i]),
                 gene_set = current_gene_set,
                 stringsAsFactors = F)
  } else {
    set_annot_df <-
      rbind(set_annot_df,
            data.frame(type = as.character(gene_set_list[i]),
                       gene_set = current_gene_set,
                       stringsAsFactors = F))
  }
  # df_filtered_gene_loc <-
  #   df_NCBI_37.3[df_NCBI_37.3$Gene_Symbol %in% current_gene_set, ]
  # write.table(df_filtered_gene_loc,
  #             file = paste("gene_loc/",
  #                          str_replace_all(gene_set_list[i],
  #                                          pattern = "GeneSet",
  #                                          replacement = "gene_loc"),
  #                          sep = ""),
  #             quote = F, sep = "\t",
  #             col.names = F, row.names = F)
}

write.table(set_annot_df,
            file = "gene_loc/annot_100k_20k/full_DE_set_annot_file.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = F)
