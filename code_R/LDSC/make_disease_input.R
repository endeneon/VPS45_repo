# 18 Mar 2022
# make the input file of each disease (.ldct)

# init
library(stringr)
library(readr)

ld_files_list <-
  read_delim("ldscore_file_list.ldcts", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
ld_files_list <- unlist(ld_files_list)

# make name column
ld_files_name <- 
  unlist(str_split(string = ld_files_list,
                   pattern = "_sig",
                   simplify = T)[, 1])
ld_control_files <-
  paste0(',ldscore/',
         unlist(str_split(string = ld_files_name,
                          pattern = "_",
                          simplify = T)[, 1]),
         '_all.')
ld_files_list_long <-
  str_c('ldscore/',
        ld_files_list,
        ld_control_files)

ld_files_df <-
  data.frame(samples = ld_files_name,
             path = ld_files_list_long)

write.table(ld_files_df,
            file = "Mar18_4_cell_types_gene_ldsc.ldcts",
            quote = F, sep = "\t",
            row.names = F, col.names = F)


ld_files_list_long2 <-
  str_c('ldscore/',
        ld_files_list)

ld_files_conc <- 
  paste(ld_files_list_long2,
        collapse = ',')
write.table(ld_files_conc,
            file = "16_types_conc",
            quote = F, sep = "\t",
            row.names = F, col.names = F)
