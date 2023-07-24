# Siwei 29 Sept 2020
# calc correlation of CROPSeq neuron

# init
library(readr)

#
VPS45_table_CROPSeq_3_lines_edgeR$Chr <- NULL
VPS45_table_CROPSeq_3_lines_edgeR$Start <- NULL
VPS45_table_CROPSeq_3_lines_edgeR$End <- NULL
VPS45_table_CROPSeq_3_lines_edgeR$Strand <- NULL

VPS45_table_CROPSeq_3_lines_edgeR$Length <- NULL

save.image(file = "VPS45_table_CROPSeq_3_lines_edgeR.RData")
