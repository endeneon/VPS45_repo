# plot peak file in VPS45 AA vs GG isogenic CD11 lines
# Siwei 04 Jun 2021
# use bigWig files as data source

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
##########
options(ucscChromosomeNames = F)
##########

# function for plotting
plot_anywhere <- function(chr, start, end, gene_name = "",
                          SNPname = "", SNPposition = 1L,
                          horizon.scale = 250,
                          mcols = 100, strand = "+",
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
  # dTrack_AA_1 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_AA_1.bigWig",
  #                        genome = "hg38", 
  #                        baseline = 0, type = "mountain", 
  #                        span = 0.1, evaluation = 1000, degree = 4, family = "gaussian",
  #                        col.title = "black",
  #                        col.axis = "black", rotation.title = 90, 
  #                        fontcolor = "black", fontsize = 16,
  #                        chromosome = chr,
  #                        start = start,
  #                        end = end,
  #                        name = "AA-1")
  dTrack_AA_1 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_AA_1.bigWig",
                           genome = "hg38", type = "horizon", 
                           col.title = "black",
                           col.axis = "black", rotation.title = 90, 
                           fontcolor = "black", fontsize = 16,
                           horizon.origin = 0, horizon.scale = horizon.scale,
                           col.horizon = "darkred",
                           fill.horizon = rev(c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3")), #rep_len("darkred", length.out = 6),
                           col.regions = c("darkred", "darkred", "darkred"),
                           chromosome = chr,
                           start = start,
                           end = end,
                           name = "AA-1")
  dTrack_AA_2 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_AA_2.bigWig",
                         genome = "hg38", type = "horizon", 
                         col.title = "black",
                         col.axis = "black", rotation.title = 90,
                         fontcolor = "black", fontsize = 16,
                         horizon.origin = 0, horizon.scale = horizon.scale,
                         col.horizon = "darkred",
                         fill.horizon = rev(c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3")), #rep_len("darkred", length.out = 6),
                         col.regions = c("darkred", "darkred", "darkred"),
                         chromosome = chr,
                         start = start,
                         end = end,
                         name = "AA-2")
  dTrack_GG_1 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_GG_1.bigWig",
                         genome = "hg38", type = "horizon", 
                         col.title = "black",
                         col.axis = "black", rotation.title = 90,
                         fontcolor = "black", fontsize = 16,
                         horizon.origin = 0, horizon.scale = horizon.scale,
                         col.horizon = "darkblue",
                         fill.horizon = c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3"), #rep_len("darkred", length.out = 6),
                         col.regions = c("darkblue", "darkblue", "darkblue"),
                         chromosome = chr,
                         start = start,
                         end = end,
                         name = "GG-1")
  dTrack_GG_2 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_GG_2.bigWig",
                         genome = "hg38", type = "horizon", 
                         col.title = "black",
                         col.axis = "black", rotation.title = 90,
                         fontcolor = "black", fontsize = 16,
                         horizon.origin = 0, horizon.scale = horizon.scale,
                         col.horizon = "darkblue",
                         fill.horizon = c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3"), #rep_len("darkred", length.out = 6),
                         col.regions = c("darkblue", "darkblue", "darkblue"),
                         chromosome = chr,
                         start = start,
                         end = end,
                         name = "GG-2")

  
  snpTrack <- AnnotationTrack(start = SNPposition, end = SNPposition, 
                              chromosome = chr,
                              id = SNPname, shape = "box",
                              name = , strand = "*",
                              group = SNPname,
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")
  
  ########### plotting
  ucscGenes <- UcscTrack(genome = genome(cell_type), table = "ncbiRefSeq",
                         track = 'NCBI RefSeq', trackType = "GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
                         rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name",
                         strand = "strand", stacking = 'pack', 
                         showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- 
    as.vector(mapIds(org.Hs.eg.db, 
                     gsub("\\.[1-9]$", "",
                          mcols(z)$symbol), "SYMBOL","REFSEQ"))
  grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor <- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 18
  
  htTrack <- HighlightTrack(trackList = list(dTrack_AA_1, dTrack_AA_2,
                                             dTrack_GG_1, dTrack_GG_2),
                            start = c(150066791, 150268163), 
                            end = c(150068025, 150271217),
                            chromosome = "chr1",
                            col = "darkorange3",
                            fill = "yellow",
                            inBackground = T)
  

  plotTracks(list(iTrack, gTrack, 
                  # htTrack,
                  dTrack_AA_1, dTrack_AA_2,
                  dTrack_GG_1, dTrack_GG_2,
                  grTrack, snpTrack),
             sizes = c(0.5, 0.5, 
                       0.5, 0.5, 0.5, 0.5,
                       0.5, 0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,
  
}


####
plot_anywhere(chr = "chr1",
              start = 150066791,
              end = 150068025,
              SNPname = "rs2027349",
              SNPposition = 150067621,
              horizon.scale = 250)

plot_anywhere(chr = "chr1",
              start = 150268163,
              end = 150275217,
              SNPname = "rs73011383",
              SNPposition = 150272921,
              horizon.scale = 250)
                  
###

# function for plotting
plot_anywhere_filled <- function(chr, start, end, gene_name = "",
                          SNPname = "", SNPposition = 1L,
                          horizon.scale = 250,
                          mcols = 100, strand = "+",
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
  # dTrack_AA_1 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_AA_1.bigWig",
  #                        genome = "hg38", 
  #                        baseline = 0, type = "mountain", 
  #                        span = 0.1, evaluation = 1000, degree = 4, family = "gaussian",
  #                        col.title = "black",
  #                        col.axis = "black", rotation.title = 90, 
  #                        fontcolor = "black", fontsize = 16,
  #                        chromosome = chr,
  #                        start = start,
  #                        end = end,
  #                        name = "AA-1")
  dTrack_AA_1 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_AA_1.bigWig",
                           genome = "hg38", type = "horizon", 
                           col.title = "black",
                           col.axis = "black", rotation.title = 90, 
                           fontcolor = "black", fontsize = 16,
                           horizon.origin = 0, horizon.scale = horizon.scale,
                           col.horizon = "darkred",
                           fill.horizon = rep_len("darkred", length.out = 6),
                           col.regions = c("darkred", "darkred", "darkred"),
                           chromosome = chr,
                           start = start,
                           end = end,
                           name = "AA-1")
  dTrack_AA_2 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_AA_2.bigWig",
                           genome = "hg38", type = "horizon", 
                           col.title = "black",
                           col.axis = "black", rotation.title = 90,
                           fontcolor = "black", fontsize = 16,
                           horizon.origin = 0, horizon.scale = horizon.scale,
                           col.horizon = "darkred",
                           fill.horizon = rep_len("darkred", length.out = 6),
                           col.regions = c("darkred", "darkred", "darkred"),
                           chromosome = chr,
                           start = start,
                           end = end,
                           name = "AA-2")
  dTrack_GG_1 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_GG_1.bigWig",
                           genome = "hg38", type = "horizon", 
                           col.title = "black",
                           col.axis = "black", rotation.title = 90,
                           fontcolor = "black", fontsize = 16,
                           horizon.origin = 0, horizon.scale = horizon.scale,
                           col.horizon = "darkblue",
                           fill.horizon = rep_len("darkblue", length.out = 6),
                           col.regions = c("darkblue", "darkblue", "darkblue"),
                           chromosome = chr,
                           start = start,
                           end = end,
                           name = "GG-2")
  dTrack_GG_2 <- DataTrack(range = "~/Data/FASTQ/Duan_Project_016/Duan_Project_016/VPS45_NPC/pre_bams/WASPed_BAMs/bigWig/R_GG_2.bigWig",
                           genome = "hg38", type = "horizon", 
                           col.title = "black",
                           col.axis = "black", rotation.title = 90,
                           fontcolor = "black", fontsize = 16,
                           horizon.origin = 0, horizon.scale = horizon.scale,
                           col.horizon = "darkblue",
                           fill.horizon = rep_len("darkblue", length.out = 6),
                           col.regions = c("darkblue", "darkblue", "darkblue"),
                           chromosome = chr,
                           start = start,
                           end = end,
                           name = "GG-2")
  
  
  snpTrack <- AnnotationTrack(start = SNPposition, end = SNPposition, 
                              chromosome = chr,
                              id = SNPname, shape = "box",
                              name = , strand = "*",
                              group = SNPname,
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")
  
  ########### plotting
  ucscGenes <- UcscTrack(genome = genome(cell_type), table = "ncbiRefSeq",
                         track = 'NCBI RefSeq', trackType = "GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
                         rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name",
                         strand = "strand", stacking = 'pack', 
                         showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- 
    as.vector(mapIds(org.Hs.eg.db, 
                     gsub("\\.[1-9]$", "",
                          mcols(z)$symbol), "SYMBOL","REFSEQ"))
  grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor <- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 18
  
  htTrack <- HighlightTrack(trackList = list(dTrack_AA_1, dTrack_AA_2,
                                             dTrack_GG_1, dTrack_GG_2),
                            start = c(150066791, 150268163), 
                            end = c(150068025, 150271217),
                            chromosome = "chr1",
                            col = "darkorange3",
                            fill = "yellow",
                            inBackground = T)
  
  plotTracks(list(iTrack, gTrack,
                  htTrack,
                  grTrack, snpTrack),
             sizes = c(0.5, 0.5,
                       0.5, 0.5, 0.5, 0.5,
                       0.5, 0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,
  
  # plotTracks(list(iTrack, gTrack,
  #                 dTrack_AA_1, dTrack_AA_2,
  #                 dTrack_GG_1, dTrack_GG_2,
  #                 grTrack, snpTrack),
  #            sizes = c(0.5, 0.5,
  #                      0.5, 0.5, 0.5, 0.5,
  #                      2, 0.5),
  #            chromosome = cell_type@ranges@NAMES,
  #            from = (cell_type@ranges@start - x_offset_1),
  #            to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
  #            transcriptAnnotation = "transcript",
  #            collapseTranscripts = "transcript")#,
  
}

####
plot_anywhere_filled(chr = "chr1",
                     start = 150040342,
                     end = 150277890,
                     SNPname = "rs2027349",
                     SNPposition = 150067621,
                     horizon.scale = 350)

plot_anywhere_filled(chr = "chr1",
                     start = 150015497,
                     end = 150504982,
                     SNPname = "rs2027349",
                     SNPposition = 150067621,
                     horizon.scale = 350)

plot_anywhere_filled(chr = "chr1",
                     start = 150013189,
                     end = 150747417,
                     SNPname = "rs2027349",
                     SNPposition = 150067621,
                     horizon.scale = 350)

plot_anywhere_filled(chr = "chr1",
                     start = 140013189,
                     end = 160747417,
                     SNPname = "rs2027349",
                     SNPposition = 150067621,
                     horizon.scale = 350)

plot_anywhere_filled(chr = "chr1",
                     start = 40013189,
                     end = 50747417,
                     SNPname = "rs2027349",
                     SNPposition = 150067621,
                     horizon.scale = 350)
