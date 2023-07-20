# 27 Jun 2021
# Gviz to plot Fig 2C

# init
library(Gviz)
library(ggplot2)

library(readr)
##
bgFile <- 
  system.file("extdata/test.bedGraph", 
              package = "Gviz")
testRange <- bgFile

###
# list all files in /Fig2C (VPS45)
import_file_list <- 
  list.files(path = "Fig2C/",
             all.files = F,
             recursive = F)

# load data
VPS45_data_list <- 
  vector(mode = "list",
         length = length(import_file_list))

proximal_gene_list <- 
  scan(file = "Fig2C/proximal_genes.txt",
       what = character(),
       quote = "")
# sort and unique
proximal_gene_list <- unique(sort(proximal_gene_list))



i <- 1
for (i in 1:length(import_file_list)) {
  VPS45_data_list[[i]] <-
    read_delim(paste("Fig2C/",
                     import_file_list[i],
                     sep = ""), 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE)
  ## select only genes in the list
  VPS45_data_list[[i]] <- 
    VPS45_data_list[[i]][VPS45_data_list[[i]]$Gene_Symbol %in% 
                           proximal_gene_list, ]
}

## construct Gviz instance
VPS45_plot <-
  DataTrack(start = VPS45_data_list[[1]]$start,
            end = VPS45_data_list[[1]]$stop,
            chromosome = "chr1",
            data = as.matrix(rbind(VPS45_data_list[[1]]$logFC,
                                   VPS45_data_list[[2]]$logFC,
                                   VPS45_data_list[[3]]$logFC)),
            genome = "hg38",
            groups = c("VPS45-1", 
                       "VPS45-2", 
                       "VPS45-3"),
            type = c("p", "a"),
            legend = T)

plotTracks(VPS45_plot)

min(VPS45_data_list[[1]]$start) #149884459
max(VPS45_data_list[[1]]$stop) #150513789
