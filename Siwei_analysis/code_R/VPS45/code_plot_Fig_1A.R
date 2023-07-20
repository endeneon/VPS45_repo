# 14 Jul 2021
# function for plotting NGN2, D15, and NPC tracs
# revised from plot_anywhere and plot_ASoC_composite
# Use OverlayTrack to combine data tracks
# VPS45 paper Fig 1A

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
###
options(ucscChromosomeNames = F)

# function for plotting
# note that currently using .bigWig as data input,
# Hence the track type should be DataTrack  
color.gradient <- 
  function(x, 
           colors=c("darkblue","red"), 
           colsteps=100) {
    return(colorRampPalette(colors) (colsteps) [findInterval(x, 
                                                             seq(min(x),
                                                                 max(x), 
                                                                 length.out = colsteps))])
  }

plot_combined_tracks <- function(chr, start, end, gene_name = "",
                                 SNPname = "", 
                                 SNPposition = 1L,
                                 mcols = 100, 
                                 strand = "+",
                                 x_offset_1 = 0, x_offset_2 = 0, ylimit = 800) {
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
                            fontsize = 16,
                            scale = 0.1)
  GWASTrack <- DataTrack(range = "PGC3_GWAS/PGC3_SCZ_w3_chr1_hg38.2021.bed",
                         genome = "hg38", type = "p",
                         chromosome = "chr1",
                         name = "PGC3 SCZ",
                         transformation = function(x) {
                           log_transformed <- 0 - log10(x)
                           return(log_transformed)
                         },
                         # col = color.gradient(function(x) {
                         #   log_transformed <- 0 - log10(x)
                         #   return(log_transformed)
                         # }),
                         col = "darkgreen",
                         background.title = "darkgreen",
                         fill.coverage = "darkgreen",
                         col.coverage = "darkgreen",
                         cex.axis = 1.2)
  alTrack_NGN2 <- DataTrack(range = "/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/Glut_NGN2_20_lines_30M_merged_sorted.bigWig",
                            chromosome = as.character(unique(seqnames(cell_type))),
                            genome = "hg38", type = "h",
                            background.title = "darkblue",
                            fill.coverage = "darkblue",
                            col.coverage = "darkblue",
                            col = "darkblue",
                            ylim = c(0, ylimit),
                            cex.axis = 1.2,
                            name = "NGN2-Glut")
  alTrack_D15 <- DataTrack(range = "/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/CN_D15_merged.bigWig",
                           chromosome = as.character(unique(seqnames(cell_type))),
                           genome = "hg38", type = "h",
                           background.title = "goldenrod",
                           fill.coverage = "goldenrod",
                           col.coverage = "goldenrod",
                           col = "goldenrod",
                           ylim = c(0, ylimit),
                           cex.axis = 1.2,
                           name = "iN-Glut")
  alTrack_NPC <- DataTrack(range = "/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/NPC_merged.bigWig",
                           chromosome = as.character(unique(seqnames(cell_type))),
                           genome = "hg38", type = "h",
                           background.title = "maroon",
                           fill.coverage = "maroon",
                           col.coverage = "maroon",
                           col = "maroon",
                           ylim = c(0, ylimit),
                           cex.axis = 1.2,
                           name = "NPC")

  
  snpTrack <- AnnotationTrack(start = SNPposition, end = SNPposition, chromosome = chr,
                              id = SNPname, shape = "box",
                              name = SNPname, strand = "*",
                              group = SNPname,
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")
  
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
                                          "SYMBOL","REFSEQ"))
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
  
  
  
  plotTracks(list(iTrack, gTrack, GWASTrack, 
                  alTrack_NGN2, alTrack_D15, alTrack_NPC, 
                  snpTrack, grTrack),
             sizes = c(0.5, 0.5, 1,
                       1, 1, 1,
                       0.5, 1.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,
  
}

dev.off()

plot_combined_tracks(chr = "chr1",
                     start = 150045621,
                     end = 150070621,
                     ylimit = 150, 
                     SNPname = "rs2027349", 
                     SNPposition = 150067621)

plot_combined_tracks(chr = "chr1",
                     start = 149990000,
                     end = 150278370,
                     ylimit = 300, 
                     SNPname = "rs2027349", 
                     SNPposition = 150067621)
