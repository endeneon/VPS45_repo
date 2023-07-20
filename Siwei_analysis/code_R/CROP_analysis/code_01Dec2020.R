# Siwei 01 Dec 2020
# Use 09+11 samples only
# count matrix regenerated using CellRanger 5.0
# extract cluster-specific cells for analysis
# do not use SCTransform()

# Use Nov30.RData
# make gRNA-specific analysis


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

# load data
load("~/NVME/CROPSeq_neuron_new_analysis_23Nov2020/R_0911_new_analysis/Nov30.RData")
# remove downstream analysis outputs
rm(list = c("Glut_neuron", "gRNA_list", "i",
            "DGE_raw", "DGE_design", "DGE_glmFit", "DGE_glmFTest"))


# assign gRNA identity by gene to each cell as meta.data
# merge neg-CTRL and neg-EGFP as one "neg_all_ctrl"
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- NULL
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <- 
  NGN2_neuron_gRNA_stripped@meta.data$gRNA.indiv
NGN2_neuron_gRNA_stripped$gRNA.bygRNA[str_detect(string = NGN2_neuron_gRNA_stripped$gRNA.bygRNA, 
                                                 pattern = "neg-")] <- "neg_all_ctrl"

## not merging gRNA by genes
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <-
#   as.vector(gsub("\\..*", "", NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA))
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <-
#   as.vector(gsub("-gene", "", NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA))
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA <-
#   as.vector(gsub("-\\d", "", NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygRNA))






# check plots
DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap", 
        label = T)
DimPlot(NGN2_neuron_gRNA_stripped, reduction = "tsne")

FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("MAP2", "SLC17A6", "SLC17A7"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SLC17A7", "SLC17A6"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SLC17A7"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SLC17A6"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SOX2", "GAD1"),
            blend = F)

table(NGN2_neuron_gRNA_stripped@meta.data$seurat_clusters)

### Glut group: # 0, 1, 2, 3, 7, 8, 9, 10, 13 (SLC17A6 + SLC17A7)
### vst transformed data are meant to be used for clustering only, 
### not for DE analysis, either for DESeq or EdgeR
Glut_neuron <-
  NGN2_neuron_gRNA_stripped[, NGN2_neuron_gRNA_stripped$seurat_clusters %in%
                              c(0, 1, 2, 3, 7, 8, 9, 10, 13)]
# save.image()
# Run EdgeR analysis
gRNA_list <- unique(Glut_neuron$gRNA.bygRNA)
gRNA_list <- gRNA_list[order(gRNA_list)]

# Make DGEList object
DGE_raw <- DGEList(counts = as.matrix(Glut_neuron@assays$RNA@counts),
                   group = as.factor(Glut_neuron$gRNA.bygRNA),
                   remove.zeros = T, 
                   genes = Glut_neuron@assays$RNA@data@Dimnames[1])
DGE_raw <- calcNormFactors(DGE_raw)
DGE_raw <- estimateDisp(DGE_raw)

DGE_design <- model.matrix(~ 0 + as.factor(Glut_neuron$gRNA.bygRNA),
                           data = DGE_raw$samples)
DGE_glmFit <- glmQLFit(DGE_raw,
                       design = DGE_design)

# write a for() loop walking through all genes
# control group always serve as the baseline
# neg-ctrl is at # 37
# total length of gene_list is 72

i <- 1L

for (i in 1:72) {
  try({
    print(paste(i, gRNA_list[i], sep = ", "))
    if (i %in% c(16, 37, 68)) {
      i <- i + 1
    } else if (i < 37) {
      DGE_glmFTest <- glmQLFTest(DGE_glmFit,
                                 contrast = c(rep_len(0, length.out = i - 1), 1, # control
                                              rep_len(0, length.out = 36 - i), -1, 
                                              rep_len(0, length.out = 35)))
      DGE_glmFTest <- DGE_glmFTest$table
      DGE_glmFTest$FDR <- p.adjust(DGE_glmFTest$PValue,
                                   method = "fdr")
      DGE_glmFTest$Gene_Symbol <- rownames(DGE_glmFTest)
      DGE_glmFTest <- join(DGE_glmFTest,
                           ENSG_coord_gene_gencodev28,
                           by = "Gene_Symbol", 
                           match = "first")
    } else if (i > 37) {
      DGE_glmFTest <- glmQLFTest(DGE_glmFit,
                                 contrast = c(rep_len(0, length.out = 36), -1, # control
                                              rep_len(0, length.out = i - 38), 1, 
                                              rep_len(0, length.out = length(gRNA_list) - i))) # 25
      DGE_glmFTest <- DGE_glmFTest$table
      DGE_glmFTest$FDR <- p.adjust(DGE_glmFTest$PValue,
                                   method = "fdr")
      DGE_glmFTest$Gene_Symbol <- rownames(DGE_glmFTest)
      DGE_glmFTest <- join(DGE_glmFTest,
                           ENSG_coord_gene_gencodev28,
                           by = "Gene_Symbol", 
                           match = "first")
    }
    write.table(DGE_glmFTest, 
                file = paste("by_gRNA/", gRNA_list[i], "_table.txt", sep = ""), 
                row.names = F,
                quote = F, sep = "\t")
  })
}

save.image()


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
FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("MAP2"),
            blend = F)
FeaturePlot(NGN2_neuron_with_gRNA, 
            features = c("SLC17A6",
                         "SLC17A7"),
            blend = T,
            ncol = 2, 
            cols = c("magenta", "blue"))

StackedVlnPlot(obj = NGN2_neuron_with_gRNA,
               features = c("DLG4", 
                            "SLC17A6", "SLC17A7", "GLS",
                            "GAD2", 
                            "TH", 
                            "VIM", "SOX2")) +
  coord_flip()


## assign gRNA identity to each cell
NGN2_neuron_gRNA_stripped <- NGN2_neuron_with_gRNA

NGN2_neuron_gRNA_stripped$gRNA.identity <- 
  rownames(NGN2_neuron_gRNA_only)[colMaxs(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))]

## plot VPS452/3 gRNA
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA <- "unassigned"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "VPS45-1-gene"] <- "VPS45_gRNA"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "VPS45-2-gene"] <- "VPS45_gRNA"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "VPS45-3-gene"] <- "VPS45_gRNA"

NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-CTRL-00018-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-CTRL-00022-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-1-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-2-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-3-gene"] <- "ctrl"

NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA <-
  as.factor(NGN2_neuron_gRNA_stripped$plot_VPS45_2_3_gRNA)

# set active identity
Idents(object = NGN2_neuron_gRNA_stripped) <- "plot_VPS45_2_3_gRNA"


DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        group.by = "plot_VPS45_2_3_gRNA",
        # cols = c("lightgrey", "blue", "red", "red"),
        order = c("unassigned", "ctrl",
                  "VPS45_gRNA"),
        cells.highlight = list(WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("ctrl")),
                               WhichCells(NGN2_neuron_gRNA_stripped,
                                          idents = c("VPS45_gRNA"))),
        # shuffle = T,
        sizes.highlight = .5,
        cols.highlight = c("red", "blue"))


DimPlot(NGN2_neuron_gRNA_stripped, 
        reduction = "umap",
        group.by = "plot_VPS45_2_3_gRNA",
        # cols = c("lightgrey", "blue", "red", "red"),
        order = c("unassigned", "ctrl",
                  "VPS45_gRNA"),
        cells.highlight = list(WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("ctrl")),
                               WhichCells(NGN2_neuron_gRNA_stripped, 
                                          idents = c("VPS45_2_gRNA")),
                               WhichCells(NGN2_neuron_gRNA_stripped,
                                          idents = c("VPS45_3_gRNA"))),
        # shuffle = T,
        sizes.highlight = .5,
        cols.highlight = c("red", "red", "blue"))

##### Fig S3
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("GAD2"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("DLG4"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("GLS"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("TH"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("VIM"),
            blend = F)
FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SOX2"),
            blend = F)

FeaturePlot(NGN2_neuron_gRNA_stripped, 
            features = c("SLC17A6"),
            blend = F)

#### plot cells assigned to unique/non-unique gRNAs
# make a data frame with each cell assigned to one unique gRNA identity
gRNA_stacked_plot_nonunique <- 
  rowSums(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))

# Select cells that have one gRNA expression more than 3x of all others
# 6658 cells passed this filter
NGN2_neuron_seurat_split <- NGN2_neuron_with_gRNA
NGN2_neuron_seurat_split <- NGN2_neuron_seurat_split[, 
                                                     (colMaxs(as.matrix(NGN2_neuron_seurat_split@assays$RNA@counts[
                                                       ((NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1] - 75):
                                                          NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1]), ]), value = T) >
                                                        Matrix::colSums(NGN2_neuron_seurat_split@assays$RNA@counts[
                                                          ((NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1] - 75):
                                                             NGN2_neuron_seurat_split@assays$RNA@counts@Dim[1]), ]) * 3 / 4)]
ncol(NGN2_neuron_seurat_split) 
NGN2_neuron_split_gRNA_only <- 
  subset(x = NGN2_neuron_seurat_split, 
         features = rownames(NGN2_neuron_seurat_split)[(nrow(NGN2_neuron_seurat_split) - 75):
                                                         nrow(NGN2_neuron_seurat_split)])

gRNA_stacked_plot_unique <- 
  rowSums(as.matrix(NGN2_neuron_split_gRNA_only@assays$RNA@counts))


######
Idents(object = NGN2_neuron_gRNA_stripped) <- "seurat_clusters"
DimPlot(NGN2_neuron_gRNA_stripped,
        reduction = "umap",
        repel = T,
        label = T)

Glut_subset <- 
  subset(NGN2_neuron_gRNA_stripped,
         idents = c("2","7","5","4","6","3"))

Glut_subset <- NGN2_neuron_gRNA_stripped
### make the df for testing

df_test <- 
  data.frame(category.name = Glut_subset$plot_VPS45_2_3_gRNA[Glut_subset$plot_VPS45_2_3_gRNA %in% c("ctrl", "VPS45_gRNA")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$plot_VPS45_2_3_gRNA %in% c("ctrl", "VPS45_gRNA")],
             stringsAsFactors = F)

transformed_read_matrix <-
  as.matrix(Glut_subset@assays$SCT@data)
transformed_read_matrix <- 
  as.data.frame(transformed_read_matrix)

df_test$VPS45.read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "VPS45", 
                                 Glut_subset$plot_VPS45_2_3_gRNA %in% c("ctrl", "VPS45_gRNA")])

df_test$gRNA.name <- as.factor(df_test$gRNA.name)

ggplot(df_test,
       aes(x = category.name,
           y = VPS45.read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic()

df_to_plot <- 
  data.frame(value = c(df_test$VPS45.read.value[df_test$category.name %in% "ctrl"],
                       df_test$VPS45.read.value[!(df_test$category.name %in% "ctrl")],
                       df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-2-gene"],
                       df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-3-gene"]),
             category = c(as.character(df_test$category.name)[df_test$category.name %in% "ctrl"],
                          as.character(df_test$category.name)[!(df_test$category.name %in% "ctrl")],
                          as.character(df_test$gRNA.name)[df_test$gRNA.name %in% "VPS45-2-gene"],
                          as.character(df_test$gRNA.name)[df_test$gRNA.name %in% "VPS45-3-gene"]))

df_to_plot <- 
  data.frame(value = c(df_test$VPS45.read.value[df_test$category.name %in% "ctrl"],
                       df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-1-gene"],
                       df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-2-gene"],
                       df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-3-gene"]),
             category = c(as.character(df_test$category.name)[df_test$category.name %in% "ctrl"],
                          as.character(df_test$gRNA.name)[df_test$gRNA.name %in% "VPS45-1-gene"],
                          as.character(df_test$gRNA.name)[df_test$gRNA.name %in% "VPS45-2-gene"],
                          as.character(df_test$gRNA.name)[df_test$gRNA.name %in% "VPS45-3-gene"]))


ggplot(df_to_plot,
      aes(x = category,
          y = value,
          fill = category,
          group = factor(category))) +
  geom_violin(trim = F, na.rm = T) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               fill = "grey30",
               dotsize = 0.2,
               binwidth = 0.02,
               position = "jitter",
               width = -1, 
               stackratio = 0.2) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  scale_fill_discrete(name = "",
                      type = brewer.pal(4, "Set1")) +
  ylab("Normalized expression value") +
  ylim(-0.5, 3.5) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,
                                   angle = 315,
                                   vjust = 0,
                                   hjust = 0))


write.table(df_to_plot,
            file = "VPS45_gRNA_stat.txt",
            quote = F, sep = "\t",
            col.names = T, row.names = F)

t.test(x = df_test$VPS45.read.value[df_test$category.name %in% "ctrl"],
       y = df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-1-gene"],
       alternative = "t", paired = F, var.equal = F)

t.test(x = df_test$VPS45.read.value[df_test$category.name %in% "ctrl"],
       y = df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-2-gene"],
       alternative = "t", paired = F, var.equal = F)

t.test(x = df_test$VPS45.read.value[df_test$category.name %in% "ctrl"],
       y = df_test$VPS45.read.value[df_test$gRNA.name %in% "VPS45-3-gene"],
       alternative = "t", paired = F, var.equal = F)


### test BAG5-2 gRNA

NGN2_neuron_gRNA_stripped$plot_BAG5_2_gRNA <- "unassigned"

NGN2_neuron_gRNA_stripped$plot_BAG5_2_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-CTRL-00018-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_BAG5_2_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-CTRL-00022-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_BAG5_2_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-1-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_BAG5_2_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "neg-EGFP-2-gene"] <- "ctrl"
NGN2_neuron_gRNA_stripped$plot_BAG5_2_gRNA[NGN2_neuron_gRNA_stripped$gRNA.identity %in% 
                                                "BAG5-2-gene"] <- "BAG5_2_gRNA"
Glut_subset <- 
  subset(NGN2_neuron_gRNA_stripped,
         idents = c("2","7","5","4","6","3"))

Glut_subset <- NGN2_neuron_gRNA_stripped

df_test <- 
  data.frame(category.name = Glut_subset$plot_BAG5_2_gRNA[Glut_subset$plot_BAG5_2_gRNA %in% c("ctrl", "BAG5_2_gRNA")],
             gRNA.name = Glut_subset$gRNA.identity[Glut_subset$plot_BAG5_2_gRNA %in% c("ctrl", "BAG5_2_gRNA")],
             stringsAsFactors = F)

transformed_read_matrix <-
  as.matrix(Glut_subset@assays$SCT@data)
transformed_read_matrix <- 
  as.data.frame(transformed_read_matrix)

df_test$BAG5.read.value <- 
  unlist(transformed_read_matrix[rownames(transformed_read_matrix) %in% "BAG5", 
                                 Glut_subset$plot_BAG5_2_gRNA %in% c("ctrl", "BAG5_2_gRNA")])

# df_test$gRNA.name <- factor(df_test$gRNA.name,
#                             levels = c("ctrl", "BAG5_2_gRNA"))
df_test$category.name <- factor(df_test$category.name,
                                levels = c("ctrl", "BAG5_2_gRNA"))

ggplot(df_test,
       aes(x = category.name,
           y = BAG5.read.value,
           fill = category.name)) +
  geom_violin(trim = F, na.rm = T) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 21,
               fill = "yellow",
               size = 2) +
  xlab("") +
  ylab("BAG5 Normalized expression value") +
  scale_fill_manual(name = "", 
                    values = c("lightblue", "orchid")) +
  geom_boxplot(notch = F,
               # notchwidth = 0.25,
               width = 0.1,
               fill = "transparent") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 10))

# p-value = 0.033
t.test(x = df_test$BAG5.read.value[df_test$category.name %in% "ctrl"], 
       y = df_test$BAG5.read.value[df_test$category.name %in% "BAG5_2_gRNA"], 
       var.equal = F,
       alternative = "t",
       paired = F)
