# 23 Jul 2021
# function for plotting gene tracks only for Fig 2C
# VPS45 paper Fig 1D

# init
library(Gviz)
# data("cpgIslands")
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(org.Hs.eg.db)
library(grDevices)
library(gridExtra)
library(GenomicRanges)

library(readr)
library(readxl)
###
options(ucscChromosomeNames = F)

# function for plotting rs2027349 site
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks
# load data

###
# list all files in /Fig2C (VPS45)
# import_file_list <- 
#   list.files(path = "Fig2C/",
#              all.files = F,
#              recursive = F)

VPS45_data_list <- 
  vector(mode = "list",
         length = 3)

proximal_gene_list <- 
  scan(file = "Fig2C/proximal_genes.txt",
       what = character(),
       quote = "")
# sort and unique
proximal_gene_list <- unique(sort(proximal_gene_list))

gene_coordination_data <-
  read_delim("Fig2C/VPS45-1-gene_table.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

gene_coordination_data <- 
  gene_coordination_data[, c(-(2:7))]

i <- 1
for (i in 1:3) {
  VPS45_data_list[[i]] <-
    read_excel("Fig2C/Fig_2C_data.xlsx", 
               sheet = i)
  VPS45_data_list[[i]] <-
    base::merge(x = VPS45_data_list[[i]],
                y = gene_coordination_data,
                by.x = "gene",
                by.y = "Gene_Symbol")
  ## select only genes in the list
  # VPS45_data_list[[i]] <- 
  #   VPS45_data_list[[i]][VPS45_data_list[[i]]$Gene_Symbol %in% 
  #                          proximal_gene_list, ]
}






plot_gene_track <- function(chr, start, end, 
                            SNP_position = 10000,
                            SNP_id = "SNP_id",
                            mcols = 100, strand = "+",
                            x_offset_1 = 0, x_offset_2 = 0, 
                            title_name = "") {
  # cell.type.selector <- 1L

  
  cell_type <-
    GRanges(seqnames = Rle(chr),
            seqinfo = Seqinfo(seqnames = chr,
                              genome = "hg38"),
            ranges = IRanges(start = start,
                             end = end,
                             names = chr),
            mcols = as.data.frame(mcols),
            strand = Rle(strand(strand)))
  
  print(as.character(unique(seqnames(cell_type))))
  
  iTrack <- IdeogramTrack(genome = genome(cell_type),
                          chromosome = as.character(unique(seqnames(cell_type))),
                          fontcolor = "black",
                          fontsize = 18)
  
  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 14,
                            scale = 0.1)
  

  snpTrack <- AnnotationTrack(start = SNP_position, end = SNP_position, 
                              chromosome = chr,
                              id = "rs2027349", shape = "box",
                              name = SNP_id, 
                              strand = "*",
                              group = c(SNP_id),
                              fontcolor.group = "black", 
                              fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")
  ### VPS45 Track
  ## construct Gviz instance
  VPS45_plot <-
    DataTrack(start = VPS45_data_list[[1]]$start,
              end = VPS45_data_list[[1]]$stop,
              chromosome = "chr1",
              data = as.matrix(rbind(VPS45_data_list[[1]]$beta,
                                     VPS45_data_list[[2]]$beta,
                                     VPS45_data_list[[3]]$beta)),
              genome = "hg38",
              groups = c("VPS45-1-gRNA", 
                         "VPS45-2-gRNA", 
                         "VPS45-3-gRNA"),
              type = c("p", "a"),
              # col = "black", 
              col.title = "black",
              col.axis = "black",
              just.group = "below",
              legend = T)
  
  
  ########### plotting
  ucscGenes <- UcscTrack(genome = genome(cell_type), 
                         table = "ncbiRefSeq",
                         track = 'NCBI RefSeq', 
                         trackType = "GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
                         rstarts = "exonStarts", 
                         rends = "exonEnds",
                         gene = "name", 
                         symbol = 'name', 
                         transcript = "name",
                         strand = "strand", 
                         stacking = 'pack', 
                         showID = T, 
                         geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, 
                                          gsub("\\.[1-9]$", 
                                               "",
                                               mcols(z)$symbol), 
                                          "SYMBOL",
                                          "REFSEQ"))
  grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor <- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 18
  
  ######
  
  # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
  #                           start = c(26212000, 26241000),
  #                           width = c(15000, 2000)
  #                           chromosome = as.character(unique(seqnames(cell_type))))
  
  
  
  plotTracks( list(iTrack, gTrack,  
                   VPS45_plot,
                   snpTrack,
                   grTrack),
              sizes = c(0.5, 0.5, 2, 0.5, 4),
              chromosome = cell_type@ranges@NAMES,
              from = (cell_type@ranges@start - x_offset_1),
              to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
              transcriptAnnotation = "transcript",
              collapseTranscripts = "transcript")#,
  
}

# plot_gene_track(chr = "chr1",
#                 start = 150067620,
#                 end = 150067621, 
#                 SNP_position = 150067621,
#                 SNP_id = "rs2027349",
#                 x_offset_1 = 250000, x_offset_2 = 250000,
#                 title_name = "NGN2-Glut")

plot_gene_track(chr = "chr1",
                start = 149850000,
                end = 150550000, 
                SNP_position = 150067621,
                SNP_id = "rs2027349",
                # x_offset_1 = 250000, x_offset_2 = 250000,
                title_name = "NGN2-Glut")

# 150067621 
min(VPS45_data_list[[1]]$start) #149884459
max(VPS45_data_list[[1]]$stop) #150513789
