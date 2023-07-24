# Siwei 28 Oct 2020
# plot Fig 3d

# init
library(readr)
library(readxl)

library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(gg.gap)

#
# Fig3d_data <- VPS45_AA_AG_GG_all_gene_table_062819_compare_with_scRNA_seq_for_Marc
# VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021_plot_only <-
#   read_delim("VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021_plot_only.txt/VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021_plot_only.txt",
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)
#
#
# Fig3d_data <- VPS45_gene_editing_AA_vs_GG_use_GG_baseline_20Sept2021_plot_only
Fig3d_data <- read_excel("VPS45_AA_AG_GG_all_gene_table_062819_compare with scRNA-seq_for Marc.xlsx",
                         sheet = 1)

Fig3d_data$logFC <- 0 - Fig3d_data$logFC.groupGG
Fig3d_data$log10PValue <- 0 - log10(Fig3d_data$PValue)
Fig3d_data$log10FDR <- 0 - log10(Fig3d_data$FDR)

Fig3d_data$Category <- "Fold Change < 1, FDR > 0.05"
Fig3d_data$Category[(Fig3d_data$logFC > 0) & (Fig3d_data$FDR > 0.05)] <-
  "Fold Change > 1, FDR > 0.05"
Fig3d_data$Category[(Fig3d_data$logFC < 0) & (Fig3d_data$FDR < 0.05)] <-
  "Fold Change < 1, FDR < 0.05"
Fig3d_data$Category[(Fig3d_data$logFC > 0) & (Fig3d_data$FDR < 0.05)] <-
  "Fold Change > 1, FDR < 0.05"
Fig3d_data$Category <- as.factor(Fig3d_data$Category)
Fig3d_data$label <- Fig3d_data$Gene_Symbol
Fig3d_data$label[!(Fig3d_data$label %in%
                   c("PLCH2", "OPCML", "GPM6A", "GRAMD1B", "GRIN2A", "MYT1L",
                     "NLGN4X", "SNAP91", "DCC", "CALN1", "NEBL", "ENOX1",
                     "CACNA1C", "SGCD",
                     "AC244033.2", "VPS45", "C1orf54"))] <- ""

Fig3d_data[Fig3d_data$label %in% "AC244033.2", ]
Fig3d_data[Fig3d_data$label %in% "C1orf54", ]
# Fig3d_data <- Fig3d_data[Fig3d_data$log10FDR < 2.5, ]
# Fig3d_data <- Fig3d_data[Fig3d_data$log10FDR > 0, ]
# Fig3d_data <- Fig3d_data[Fig3d_data$log10FC > -2, ]
# Fig3d_data <- Fig3d_data[Fig3d_data$log10FC < 2, ]


ggplot(Fig3d_data, aes(x = logFC,
                       y = log10PValue,
                       colour = Category,
                       label = label)) +
  geom_point(size = 0.2) +
  # scale_y_continuous(#data = Fig3d_data[Fig3d_data$Geneid != "", ],
  #                    trans = squash_axis(from = 0,
  #                                        to = 2,
  #                                        factor = 10)) +
  geom_point(data = Fig3d_data[Fig3d_data$label != "", ],
             colour = "black") +
  geom_hline(yintercept = 0 - log10(0.05),
             linetype = 2) +
  # geom_vline(xintercept = 0,
  #            linetype = 2) +

  xlab("log2FC") +
  xlim(-2, 2) +
  ylim(0, 5.5) +
  # ylim(0, 2.5) +
  geom_text_repel(data = Fig3d_data[Fig3d_data$label != "", ],
                  aes(label = label),
                  segment.colour = "black",

                  colour = "black",
                  box.padding = .2,
                  # segment.curvature = 1,
                  # segment.ncp = 3,
                  # segment.angle = 160,
                  min.segment.length = 0,
                  seed = 42) +

  guides(colour = guide_legend(override.aes = aes(color = NA))) +
  scale_colour_manual(values = c("#1F78B4", "#A6CEE3",
                                 "#E31A1C", "#FB9A99")) +
  theme_classic() +
  # scale_y_continuous(trans = "log2") +
  theme(legend.position = "none")


############
squash_axis <- function(from, to, factor) {
  # if(!is.na(from) & !is.na(to))
  # A transformation function that squashes the range of [from, to] by factor on a given axis

  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  #
  # Returns:
  #   A transformation called "squash_axis", which is capsulated by trans_new() function

  trans <- function(x) {
    if(any(is.na(x))) {
      return(x)
    }
    # print(paste(from, to, factor))
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to

    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    print(isq, ito)
    return(x)
  }

  inv <- function(x) {

    if(any(is.na(x))) {
      return(x)
    }
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor

    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))

    return(x)
  }

  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}







# trans <- function(x) {
#     ifelse(x > 1, x - 1, ifelse(x < -1, x + 1, x/4))
#   }
# inv <- function(x) {
#     ifelse(x > 0, x + 1, ifelse(x < 0, x - 1, x*4))
#   }
# my_trans <- trans_new("my_trans", trans, inv)


####

dat <- data.frame(group=rep(c('A', 'B', 'C', 'D'), each = 10),
                  value=c(rnorm(10), rnorm(10)+100))
ggplot(dat,aes(x=group,y=value))+
  geom_point()+
  scale_y_continuous(trans = squash_axis(5, 95, 10))
