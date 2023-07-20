# Siwei 12 Jul 2021
# make circos plot for iN-NGN2 and D15 Glut neurons credible SNPs

# init
library(circlize)
library(readxl)

library(RColorBrewer)

library(stringr)

# load data
df.NGN2.SNP <- read_excel("circo_plot_data.xlsx",
                          sheet = "NGN2")
df.D15.SNP <- read_excel("circo_plot_data.xlsx",
                         sheet = "D15")
df.Gene.list <- read_excel("circo_plot_data.xlsx",
                           sheet = "Genes")

# initiate plot, only use autosomes
circos.initializeWithIdeogram(species='hg38',
                              chromosome.index = paste0("chr", c(1:22)))

# test to generate random bed formats
bed <- generateRandomBed(nr = 50,
                        fun = function(k) sample(letters,
                                                 k,
                                                 replace = TRUE))
## convert df.Gene.list to bed format
df.BED.Gene.list <- data.frame(chr = df.Gene.list$CHR,
                               start = df.Gene.list$POS,
                               end = df.Gene.list$POS + 1,
                               Gene.name = df.Gene.list$Gene,
                               Type = df.Gene.list$Type,
                               stringsAsFactors = F)
df.BED.Gene.list$Type <- factor(df.BED.Gene.list$Type)

# convert df.NGN.SNP to region format
df.BED.NGN2.SNP <- data.frame(start = df.NGN2.SNP$POS,
                              end = df.NGN2.SNP$POS)


dev.off()
circos.initializeWithIdeogram(species='hg38',
                              chromosome.index = paste0("chr", c(1:22)),
                              plotType = NULL)
circos.genomicLabels(df.BED.Gene.list,
                     labels.column = 4,
                     side = "outside",
                     col = as.numeric(df.BED.Gene.list$Type))
circos.genomicIdeogram()
circos.genomicPoints(region = df.BED.NGN2.SNP,
                     sector.index = as.numeric(factor(df.NGN2.SNP$CHR)),
                     value = df.NGN2.SNP$`0-logFDR`)


circos.clear()

test_plot <- circos.initializeWithIdeogram(species='hg38',
                                           chromosome.index = paste0("chr", c(1:22)),
                                           plotType = NULL)
test_plot <- circos.genomicLabels(df.BED.Gene.list,
                                  labels.column = 4,
                                  side = "outside",
                                  col = as.numeric(df.BED.Gene.list$Type))
test_plot <- circos.genomicIdeogram()
circos.clear()







circos.track(ylim = c(0, 1), panel.fun = function(x, y) circos.genomicAxis())
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "bottom", direction = "inside")
})


####

dev.off()

set.seed(123)
circos.initializeWithIdeogram(species='hg38',
                              chromosome.index = paste0("chr", c(1:22)),
                              plotType = NULL)
circos.genomicLabels(df.BED.Gene.list,
                     labels.column = 4,
                     side = "outside",
                     col = as.numeric(df.BED.Gene.list$Type))
circos.track(ylim = c(0, 1),
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               circos.rect(xlim[1], 0,
                           xlim[2], 1,
                           col = rand_color(n = 1,
                                            luminosity = "light"))
               circos.text(mean(xlim),
                           mean(ylim),
                           str_replace_all(string = chr,
                                           pattern = "chr",
                                           replacement = ""),
                           cex = 0.7,
                           col = "black",
                           facing = "inside", niceFacing = TRUE)
             },
             track.height = 0.15,
             bg.border = NA)

set_track_gap(mm_h(0))
# plot NGN2 SNP list
circos.track(ylim = c(0, 0.5),
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               circos.rect(xlim[1], 0,
                           xlim[2], 0.5,
                           col = "lightblue",
                           border = NA)
               circos.points(x = df.NGN2.SNP[df.NGN2.SNP$CHR %in% chr, ]$POS,
                             y = df.NGN2.SNP[df.NGN2.SNP$CHR %in% chr, ]$`0-logFDR` /
                               max(df.NGN2.SNP$`0-logFDR`) / 2,
                             pch = 20)
             },
             track.height = 0.15,
             bg.border = NA)
set_track_gap(mm_h(0))
# plot D15 SNP list
circos.track(ylim = c(0, 0.5),
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               circos.rect(xlim[1], 0,
                           xlim[2], 0.5,
                           col = "khaki",
                           border = NA)
               circos.points(x = df.D15.SNP[df.D15.SNP$CHR %in% chr, ]$POS,
                             y = df.D15.SNP[df.D15.SNP$CHR %in% chr, ]$`0-logFDR` /
                               max(df.NGN2.SNP$`0-logFDR`) / 2,
                             pch = 20)
             },
             track.height = 0.15,
             bg.border = NA)






circos.genomicPoints(region = df.BED.NGN2.SNP,
                     panel.fun = function(x, y) {
                       print(CELL_META$sector.index)
                     },
                     value = df.NGN2.SNP$`0-logFDR`)
                     # sector.index = as.numeric(factor(df.NGN2.SNP$CHR)),
                     # value = df.NGN2.SNP$`0-logFDR`)