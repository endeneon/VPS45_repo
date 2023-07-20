# Siwei 31 Aug 2021
# Re-plot Fig 6F

# init
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(colorspace)

library(readxl)
library(stringr)

# load data
df_raw <-
  read_excel("Data_Fig6F.xlsx",
             sheet = 1)
record_Gene_Symbol <- df_raw$Gene_Symbol
df_raw$Gene_Symbol <- NULL
rownames(df_raw) <- record_Gene_Symbol

# make the plot (neuron differentiation)
df_to_plot <- as.matrix(df_raw[1:77, ])
rownames(df_to_plot) <- record_Gene_Symbol[1:77]

# heatmap.2(df_to_plot,
#           density.info = "none",
#           trace = "none",
#           col = diverge_hsv(100),
#           margins = c(8, 4))
dev.off()
heatmap.2(as.matrix(t(df_to_plot)),
          density.info = "none",
          trace = "none",
          Colv = T,
          Rowv = T,
          col = diverge_hsv(100),
          margins = c(2, 2),
          srtCol = 45,
          offsetRow = 0,
          keysize = 1,
          cexRow = 1)

dev.off()


# make the plot (synaptic)
df_to_plot <- as.matrix(df_raw[78:nrow(df_raw), ])
rownames(df_to_plot) <- record_Gene_Symbol[78:length(record_Gene_Symbol)]

# heatmap.2(df_to_plot,
#           density.info = "none",
#           trace = "none",
#           col = diverge_hsv(100),
#           margins = c(8, 4))
dev.off()
heatmap.2(as.matrix(t(df_to_plot)),
          density.info = "none",
          trace = "none",
          Colv = T,
          Rowv = F,
          col = diverge_hsv(100),
          margins = c(4, 1),
          srtCol = 45,
          offsetRow = -70,
          # keysize = 1,
          cexRow = 1)

dev.off()

# make the plot (neuron differentiation)
df_to_plot <- as.matrix(df_raw[1:30, ])
rownames(df_to_plot) <- record_Gene_Symbol[1:30]

# make the plot (synaptic)
df_to_plot <- as.matrix(df_raw[78:nrow(df_raw), ])
rownames(df_to_plot) <- record_Gene_Symbol[78:length(record_Gene_Symbol)]