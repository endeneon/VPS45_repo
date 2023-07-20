# 15 Jul 2021
# function for plotting NGN2, D15, and NPC tracs
# revised from plot_anywhere and plot_ASoC_composite
# Use OverlayTrack to combine data tracks
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
###
options(ucscChromosomeNames = F)

# function for plotting rs2027349 site
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks

## write a cell type selector
# cell.type.selector <- 1L
#
# if (cell.type.selector == 1L) {
#   print("NGN2-Glut")
#   background.bam <- "NGN2_Glut_rs2027349_het_merged_chr1.bam"
#   A.bam <- "output/NGN2.ref.A.bam"
#   G.bam <- "output/NGN2.alt.G.bam"
# } else if (cell.type.selector == 2L) {
#   print("D15-Glut")
#   background.bam <- "CN15_rs2027349_het_merged_chr1.bam"
#   A.bam <- "output/D15.ref.A.bam"
#   G.bam <- "output/D15.alt.G.bam"
# } else if (cell.type.selector == 3L) {
#   print("NPC")
#   background.bam <- "NPC_rs2027349_het_merged_chr1.bam"
#   A.bam <- "output/NPC.ref.A.bam"
#   G.bam <- "output/NPC.alt.G.bam"
# }


plot_AsoC_peaks <- function(chr, start, end, gene_name = "",
                            mcols = 100, strand = "+",
                            x_offset_1 = 0, x_offset_2 = 0, ylimit = 400,
                            cell.type.selector = 1L,
                            title_name = "") {
  # cell.type.selector <- 1L

  if (cell.type.selector == 1L) {
    print("NGN2-Glut")
    background.bam <- "NGN2_Glut_rs2027349_het_merged_chr1.bam"
    A.bam <- "output/NGN2.ref.A.bam"
    G.bam <- "output/NGN2.alt.G.bam"
  } else if (cell.type.selector == 2L) {
    print("D15-Glut")
    background.bam <- "CN15_rs2027349_het_merged_chr1.bam"
    A.bam <- "output/D15.ref.A.bam"
    G.bam <- "output/D15.alt.G.bam"
  } else if (cell.type.selector == 3L) {
    print("NPC")
    background.bam <- "NPC_rs2027349_het_merged_chr1.bam"
    A.bam <- "output/NPC.ref.A.bam"
    G.bam <- "output/NPC.alt.G.bam"
  }

  print(paste("/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/chr1_bams/",
              background.bam,
              collapse = "",
              sep = ""))
  print(paste("/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/chr1_bams/",
              G.bam,
              collapse = "", sep = ""))

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

  alTrack_Region <-
    AlignmentsTrack(range = paste("/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/chr1_bams/",
                                  background.bam,
                                  collapse = "",
                                  sep = ""),
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "palegreen3",
                    fill.coverage = "palegreen3",
                    col.coverage = "palegreen3",
                    ylim = c(0, ylimit),
                    cex.axis = 1,
                    col.title = "black",
                    col.axis = "black",
                    name = title_name,
                    alpha = 0.8,
                    fontsize = 14,
                    show.title = FALSE)

  alTrack_A <- AlignmentsTrack(paste("/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/chr1_bams/",
                                     A.bam,
                                     collapse = "", sep = ""),
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkred",
                               fill.coverage = "darkred",
                               col.coverage = "darkred",
                               ylim = c(0, ylimit),
                               cex.axis = 1,
                               # col.title="black",
                               name = "Allele_A")

  alTrack_G <- AlignmentsTrack(paste("/home/zhangs3/Data/FASTQ/Duan_Project_014/bigWig/chr1_bams/",
                                     G.bam,
                                     collapse = "", sep = ""),
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkblue",
                               fill.coverage = "darkblue",
                               col.coverage = "darkblue",
                               ylim = c(0, ylimit),
                               cex.axis = 1,
                               # col.title="black",
                               name = "Allele_G")

  OvlTrack <- OverlayTrack(trackList = list(alTrack_Region, alTrack_G, alTrack_A),
                           name = "REF_ALT_REGION",
                           fontsize = 14)

  snpTrack <- AnnotationTrack(start = 150067621, end = 150067621, chromosome = "chr1",
                              id = "rs2027349", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs2027349"),
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



  plotTracks( list(gTrack, OvlTrack,
                   snpTrack),
              sizes = c(1, 4, 0.5),
              chromosome = cell_type@ranges@NAMES,
              from = (cell_type@ranges@start - x_offset_1),
              to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
              transcriptAnnotation = "transcript",
              collapseTranscripts = "transcript")#,

}

plot_AsoC_peaks(chr = "chr1",
                start = 150067551,
                end = 150067700,
                ylimit = 700,
                cell.type.selector = 1L,
                title_name = "NGN2-Glut")

plot_AsoC_peaks(chr = "chr1",
                start = 150067551,
                end = 150067700,
                ylimit = 700,
                cell.type.selector = 2L,
                title_name = "iN-Glut")

plot_AsoC_peaks(chr = "chr1",
                start = 150067551,
                end = 150067700,
                ylimit = 700,
                cell.type.selector = 3L,
                title_name = "NPC")
