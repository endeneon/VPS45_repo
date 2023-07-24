# make additional plots for CSC revision


# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
plan("multisession", workers = 12)
# plan("sequential")
# plan()
options(expressions = 500000)
options(future.globals.maxSize = 2097152000)
set.seed(42)
### There might be conflicts between packages,
### if running SCTransform, only load the three libraries above
### SCTransform is not compatible with Rfast
##############################################################


# library(data.table)
library(ggplot2)
library(RColorBrewer)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)

library(edgeR)

library(readr)

# load data
load("Nov30.RData")
# remove downstream analysis outputs
rm(list = c("Glut_neuron", "gRNA_list", "i",
            "DGE_raw", "DGE_design", "DGE_glmFit", "DGE_glmFTest"))



DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        label = T, 
        repel = T)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("MAP2"),
            blend = F)

### plotting
NGN2_neuron_with_gRNA <- SCTransform(NGN2_neuron_gRNA_stripped,
                                     method = "glmGamPoi",
                                     vars.to.regress = "percent.mt",
                                     variable.features.n = 10000, 
                                     verbose = T)
NGN2_neuron_with_gRNA <- RunTSNE(NGN2_neuron_with_gRNA, 
                                 tsne.method = "Rtsne", 
                                 seed.use = 42)
NGN2_neuron_with_gRNA <- RunUMAP(NGN2_neuron_with_gRNA, 
                                 reduction = "pca",
                                 dims = 1:30,
                                 verbose = T, 
                                 seed.use = 42)
NGN2_neuron_with_gRNA <- FindNeighbors(NGN2_neuron_with_gRNA,
                                       k.param = 5,
                                       nn.method = "rann",
                                       reduction = "pca",
                                       dims = 1:30)
## here FindCluster() requires single-thread
NGN2_neuron_with_gRNA <- FindClusters(NGN2_neuron_with_gRNA,
                                      resolution = 0.18, 
                                      graph.name = "RNA_snn",
                                      verbose = T)
DimPlot(NGN2_neuron_with_gRNA, 
        reduction = "umap",
        label = T, 
        repel = T)

######
Idents(object = NGN2_neuron_with_gRNA) <- "seurat_clusters"
DimPlot(NGN2_neuron_with_gRNA,
        reduction = "umap",
        repel = T,
        label = T)

Glut_subset <- 
  subset(NGN2_neuron_with_gRNA,
         idents = c("2","7","5","4","6","3"))

ncol(Glut_subset)

## Add line-specific identitites
CD_09_barcodes <- 
  read_csv("CD_09_barcodes.txt", 
           col_names = FALSE)
CD_09_barcodes <-
  CD_09_barcodes$X1

CD_11_barcodes <- 
  read_csv("CD_11_barcodes.txt", 
           col_names = FALSE)
CD_11_barcodes <-
  CD_11_barcodes$X1
length(CD_11_barcodes)

Glut_subset@meta.data$line_ident <- "Undetermined"
Glut_subset@meta.data$line_ident[colnames(Glut_subset) %in% CD_09_barcodes] <- "CD-09"
Glut_subset@meta.data$line_ident[colnames(Glut_subset) %in% CD_11_barcodes] <- "CD-11"

## subset Glut_subset by cell line to plot separately #####

## line 09 #####

plot_Glut_subset <-
  Glut_subset[, Glut_subset$line_ident == "CD-09"]

transformed_read_matrix <-
  as.data.frame(as.matrix(plot_Glut_subset@assays$SCT@data))

# plot VPS45
plot_Glut_subset$plot_identity <- "unassigned"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "VPS45-2-gene"] <- "VPS45_gRNA"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "VPS45-3-gene"] <- "VPS45_gRNA"

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                                "neg-CTRL-00018-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                                "neg-CTRL-00022-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                                "neg-EGFP-1-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                                "neg-EGFP-2-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                                "neg-EGFP-3-gene"] <- "ctrl"

# plot_Glut_subset$plot_identity <- as.factor(plot_Glut_subset$plot_identity)

# Set active identity
Idents(object = plot_Glut_subset) <- "plot_identity"

# VPS45
df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "VPS45",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")]),
             stringsAsFactors = F)
# df_test$gRNA.name <- as.factor(df_test$gRNA.name)

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09, VPS45 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 7.29e-5",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "VPS45_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "VPS45_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5 1/2/3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-1-gene"] <- "BAG5_1_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                            levels = c("ctrl", "BAG5_1_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.437",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_1_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_1_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])
# BAG5-2

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-2-gene"] <- "BAG5_2_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_2_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.818",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_2_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_2_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5-3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-3-gene"] <- "BAG5_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.433",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])


### PBML1/GNL3 gRNA_3
# PBML1

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "PBRM1-3-gene"] <- "PBRM1_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "PBRM1",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "PBRM1_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09, PBRM1 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.729",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])


# GNL3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "PBRM1-3-gene"] <- "PBRM1_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "GNL3",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "PBRM1_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09, GNL3 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.468",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

save.image(file = "NGN2_plot_27Feb2023")

## line 11 #####

plot_Glut_subset <-
  Glut_subset[, Glut_subset$line_ident == "CD-11"]

transformed_read_matrix <-
  as.data.frame(as.matrix(plot_Glut_subset@assays$SCT@data))

# plot VPS45
plot_Glut_subset$plot_identity <- "unassigned"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "VPS45-2-gene"] <- "VPS45_gRNA"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "VPS45-3-gene"] <- "VPS45_gRNA"

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-CTRL-00018-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-CTRL-00022-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-EGFP-1-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-EGFP-2-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-EGFP-3-gene"] <- "ctrl"

# plot_Glut_subset$plot_identity <- as.factor(plot_Glut_subset$plot_identity)

# Set active identity
Idents(object = plot_Glut_subset) <- "plot_identity"

# VPS45
df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "VPS45",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")]),
             stringsAsFactors = F)
# df_test$gRNA.name <- as.factor(df_test$gRNA.name)

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, VPS45 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 4.44e-8",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "VPS45_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "VPS45_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5 1/2/3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-1-gene"] <- "BAG5_1_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_1_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.437",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_1_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_1_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5-2

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-2-gene"] <- "BAG5_2_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_2_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.181",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_2_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_2_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5-3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-3-gene"] <- "BAG5_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.489",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

### PBML1/GNL3 gRNA_3
# PBML1

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "PBRM1-3-gene"] <- "PBRM1_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "PBRM1",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "PBRM1_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, PBRM1 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.032",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# GNL3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "PBRM1-3-gene"] <- "PBRM1_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "GNL3",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "PBRM1_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, GNL3 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.023",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

save.image(file = "NGN2_plot_27Feb2023")


## line 09+11 #####

plot_Glut_subset <-
  Glut_subset[, Glut_subset$line_ident %in% c("CD-09", "CD-11")]

transformed_read_matrix <-
  as.data.frame(as.matrix(plot_Glut_subset@assays$SCT@data))

# plot VPS45
plot_Glut_subset$plot_identity <- "unassigned"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "VPS45-2-gene"] <- "VPS45_gRNA"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "VPS45-3-gene"] <- "VPS45_gRNA"

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-CTRL-00018-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-CTRL-00022-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-EGFP-1-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-EGFP-2-gene"] <- "ctrl"
plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "neg-EGFP-3-gene"] <- "ctrl"

# plot_Glut_subset$plot_identity <- as.factor(plot_Glut_subset$plot_identity)

# Set active identity
Idents(object = plot_Glut_subset) <- "plot_identity"

# VPS45
df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "VPS45",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "VPS45_gRNA")]),
             stringsAsFactors = F)
# df_test$gRNA.name <- as.factor(df_test$gRNA.name)

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09+11, VPS45 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 1.39e-11",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "VPS45_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "VPS45_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5 1/2/3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-1-gene"] <- "BAG5_1_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_1_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_1_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.437",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_1_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_1_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5-2

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-2-gene"] <- "BAG5_2_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_2_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_2_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.181",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_2_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_2_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# BAG5-3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "BAG5-3-gene"] <- "BAG5_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "BAG5_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-11, BAG5 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.489",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "BAG5_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "BAG5_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

### PBML1/GNL3 gRNA_3
# PBML1

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "PBRM1-3-gene"] <- "PBRM1_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "PBRM1",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "PBRM1_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09+11, PBRM1 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.038",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# GNL3

plot_Glut_subset$plot_identity[plot_Glut_subset$gRNA.indiv %in%
                                 "PBRM1-3-gene"] <- "PBRM1_3_gRNA"

df_test <- 
  data.frame(category.name = plot_Glut_subset$plot_identity[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             gRNA.name = plot_Glut_subset$gRNA.indiv[plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")],
             read.value = unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "GNL3",
                                                         plot_Glut_subset$plot_identity %in% c("ctrl", "PBRM1_3_gRNA")]),
             stringsAsFactors = F)
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "PBRM1_3_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  geom_jitter(size = 0.5,
              width = 0.1) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_manual(values = c("lightskyblue", "orchid")) +
  ylab("Normalized expression value") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12)) +
  ggtitle(paste("CD-09+11, GNL3 expression;",
                "Student's t-test (two-sided, non-parametric)",
                "P value = 0.021",
                sep = "\n"))

t.test(df_test$read.value[df_test$category.name == "ctrl"],
       df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"],
       var.equal = F, alternative = "t", paired = F)

mean(df_test$read.value[df_test$category.name == "PBRM1_3_gRNA"]) /
  mean(df_test$read.value[df_test$category.name == "ctrl"])

# save.image(file = "NGN2_plot_27Feb2023")
