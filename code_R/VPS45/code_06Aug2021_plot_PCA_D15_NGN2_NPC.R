# Siwei 06 Aug 2021
# make PCA plot for D15, NGN2, and NPCs

# init
library(readr)

library(ggplot2)
library(gplots)

library(factoextra)
library(stringr)

# load data
df_raw_input <- 
  read_delim("NGN2_D15_NPC_summit_500bp_4_PCA.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE, skip = 1)

df_read_count <- df_raw_input[, 7:ncol(df_raw_input)]

# trim excessive words in sample names (colnames)
sample_names <- colnames(df_read_count)
sample_names <- 
  str_remove_all(sample_names,
                 pattern = '_WASPed.bam')
sample_names <-
  str_remove_all(sample_names,
                 pattern = '_rapid_neuron20')
sample_names <-
  str_remove_all(sample_names,
                 pattern = "_S._new$")
sample_names <-
  str_remove_all(sample_names,
                 pattern = "_S.._new$")
# reassign colnames
colnames(df_read_count) <- sample_names

# calculate PCA use prcomp
p <- prcomp(as.matrix(t(df_read_count)),
            center = T,
            scale. = T,
            tol = 0)
scr_p <- fviz_eig(p)
fviz_pca_ind(p,
             habillage = as.factor(c(rep_len("D15", length.out = 20),
                                     rep_len("NGN2", length.out = 20),
                                     rep_len("NPC", length.out = 20))),
             palette = "Dark2",
             # label = "none",
             addEllipses = T,
             invisible = "quali",
             repel = T)


fviz_pca_var(p,
             geom.var = "text",
             label = "none")
