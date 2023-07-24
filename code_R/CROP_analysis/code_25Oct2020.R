# Siwei 25 Oct 2020
# Use 09+11 samples only
# count matrix regenerated using CellRanger 5.0
# extract cluster-specific cells for analysis
# do not use SCTransform()

# init
library(Seurat)
library(data.table)
# library(ggplot2)
# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)
library(future)


# set threads
# setDTthreads(threads = 10)
# this one will cause Cstack error if used with SCTransform() or NormalizeData()
# plan("multisession", workers = 2)
plan("sequential")

# load data use read10x_h5
NGN2_neuron_0911 <- Read10X_h5(filename = "../neuron_0911_new_23Nov2020/outs/filtered_feature_bc_matrix.h5")
NGN2_neuron_0911 <- CreateSeuratObject(counts = NGN2_neuron_0911, 
                                       project = "NGN2_neuron_0911")

nrow(NGN2_neuron_0911) # 29922 genes
ncol(NGN2_neuron_0911) # 13270 cells

# make a seurat object for independent operation
NGN2_neuron_seurat <- NGN2_neuron_0911
# Run normalisation first to separate clusters
# including NormalizeData(), ScaleData(), and FindVariableFeatures

# data preprocessing

# remove genes with all 0 count, 29922 -> 24766
NGN2_neuron_seurat <- subset(x = NGN2_neuron_seurat, 
                             features = rownames(NGN2_neuron_seurat)[rowMaxs(as.matrix(NGN2_neuron_seurat@assays$RNA@counts), value = T) > 0])
nrow(NGN2_neuron_seurat)

## remove cells without any gRNA reads
## 10247 cells with valid gRNAs
NGN2_neuron_seurat <- NGN2_neuron_seurat[, 
                                         Matrix::colSums(NGN2_neuron_seurat[(nrow(NGN2_neuron_seurat) - 75)
                                                                            :nrow(NGN2_neuron_seurat), ]) > 0]
ncol(NGN2_neuron_seurat)

# Select cells that have one gRNA expression more than 3x of all others
# 6658 cells passed this filter
NGN2_neuron_seurat <- NGN2_neuron_seurat[, 
                                         (colMaxs(as.matrix(NGN2_neuron_seurat@assays$RNA@counts[
                                           ((NGN2_neuron_seurat@assays$RNA@counts@Dim[1] - 75):
                                              NGN2_neuron_seurat@assays$RNA@counts@Dim[1]), ]), value = T) >
                                            Matrix::colSums(NGN2_neuron_seurat@assays$RNA@counts[
                                              ((NGN2_neuron_seurat@assays$RNA@counts@Dim[1] - 75):
                                                 NGN2_neuron_seurat@assays$RNA@counts@Dim[1]), ]) * 3 / 4)]
ncol(NGN2_neuron_seurat)


# Run standard workflow with NormalizeData() and FindVariableFeatures()

NGN2_neuron_seurat <- NormalizeData(NGN2_neuron_seurat,
                                    # normalization.method = "LogNormalize",
                                    # scale.factor = 10000,
                                    verbose = T)

NGN2_neuron_seurat <- PercentageFeatureSet(NGN2_neuron_seurat, 
                                           pattern = "^MT-", 
                                           col.name = "percent.mt")

NGN2_neuron_seurat <- FindVariableFeatures(NGN2_neuron_seurat, 
                                           selection.method = "vst",
                                           verbose = T)
NGN2_neuron_seurat <- ScaleData(NGN2_neuron_seurat, 
                                features = rownames(NGN2_neuron_seurat),
                                vars.to.regress = "percent.mt", 
                                verbose = T, 
                                model.use = "negbinom")
save.image()

# Run PCA, tSNE and UMAP
NGN2_neuron_seurat <- RunPCA(NGN2_neuron_seurat, 
                             verbose = T)
ElbowPlot(NGN2_neuron_seurat)
NGN2_neuron_seurat <- RunTSNE(NGN2_neuron_seurat, 
                                 dims = 1:20,
                                 verbose = T)
NGN2_neuron_seurat <- RunUMAP(NGN2_neuron_seurat, 
                              dims = 1:20,
                              verbose = T)
FeaturePlot(NGN2_neuron_seurat, 
            features = c("VPS45", "VPS45-2-gene"),
            blend = F)
#######

# split gRNA rows from Seurat matrix for subsequent quantification
NGN2_neuron_gRNA_stripped <- 
  subset(x = NGN2_neuron_seurat, 
         features = rownames(NGN2_neuron_seurat)[1:(nrow(NGN2_neuron_seurat) - 76)])

# set another Seurat matrix contains gRNA counts only
NGN2_neuron_gRNA_only <- 
  subset(x = NGN2_neuron_seurat, 
         features = rownames(NGN2_neuron_seurat)[(nrow(NGN2_neuron_seurat) - 75):nrow(NGN2_neuron_seurat)])

# make a data frame that each cell assigned to one unique gRNA identity, 6,658 cells
cell_gRNA_identity <- rownames(NGN2_neuron_gRNA_only)[colMaxs(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))]
cell_gRNA_identity <- as.data.frame(rbind(cell_gRNA_identity), 
                                    stringsAsFactors = F)
colnames(cell_gRNA_identity) <- colnames(NGN2_neuron_gRNA_only)
# check number of cells assigned to each gRNA
table(rownames(NGN2_neuron_gRNA_only)[colMaxs(as.matrix(NGN2_neuron_gRNA_only@assays$RNA@counts))])

# perform dimensionality reduction, use PCA first then t-SNE, umap
NGN2_neuron_gRNA_stripped <- FindVariableFeatures(NGN2_neuron_gRNA_stripped, 
                                                  verbose = T)
NGN2_neuron_gRNA_stripped <- RunPCA(NGN2_neuron_gRNA_stripped, 
                                    # features = 3000,
                                    verbose = T)
# check the dimensionality of signal, use PC 1-20
ElbowPlot(NGN2_neuron_gRNA_stripped)
# tsne
NGN2_neuron_gRNA_stripped <- RunTSNE(NGN2_neuron_gRNA_stripped, 
                                     dims = 1:20, 
                                     verbose = T)
NGN2_neuron_gRNA_stripped <- RunUMAP(NGN2_neuron_gRNA_stripped, 
                                     dims = 1:30,
                                     verbose = T)
# find neibours and cluster
NGN2_neuron_gRNA_stripped <- FindNeighbors(NGN2_neuron_gRNA_stripped,
                                           dims = 1:20,
                                           do.plot = T, 
                                           verbose = T)
NGN2_neuron_gRNA_stripped <- FindClusters(NGN2_neuron_gRNA_stripped, 
                                          resolution = 0.5)
# assign gRNA identity to each cell as meta.data
NGN2_neuron_gRNA_stripped@meta.data$gRNA.indiv <- 
  as.character(as.vector(cell_gRNA_identity[1, ]))

# assign gRNA identity by gene to each cell as meta.data
# merge neg-CTRL and neg-EGFP as one "neg_all_ctrl"
# NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene <- NULL
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene <- 
  NGN2_neuron_gRNA_stripped@meta.data$gRNA.indiv
NGN2_neuron_gRNA_stripped$gRNA.bygene[str_detect(string = NGN2_neuron_gRNA_stripped$gRNA.bygene, 
                                                 pattern = "neg-")] <- "neg_all_ctrl"
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene <-
  as.vector(gsub("\\..*", "", NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene))
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene <-
  as.vector(gsub("-gene", "", NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene))
NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene <-
  as.vector(gsub("-\\d", "", NGN2_neuron_gRNA_stripped@meta.data$gRNA.bygene))


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
save.image()
# Run EdgeR analysis
gRNA_list <- unique(Glut_neuron$gRNA.bygene)
gRNA_list <- gRNA_list[order(gRNA_list)]

# Make DGEList object
DGE_raw <- DGEList(counts = as.matrix(Glut_neuron@assays$RNA@counts),
                   group = as.factor(Glut_neuron$gRNA.bygene),
                   remove.zeros = T, 
                   genes = Glut_neuron@assays$RNA@data@Dimnames[1])
DGE_raw <- calcNormFactors(DGE_raw)
DGE_raw <- estimateDisp(DGE_raw)

DGE_design <- model.matrix(~ 0 + as.factor(Glut_neuron$gRNA.bygene),
                           data = DGE_raw$samples)
DGE_glmFit <- glmQLFit(DGE_raw,
                       design = DGE_design)

# write a for() loop walking through all genes
# control group always serve as the baseline
# neg-ctrl is at #13
# total length of gene_list is 25

i <- 1L

for (i in 1:25) {
  try({
    print(paste(i, gRNA_list[i], sep = ", "))
    if (i == 13) {
      i <- i + 1
    } else if (i < 13) {
      DGE_glmFTest <- glmQLFTest(DGE_glmFit,
                                 contrast = c(rep_len(0, length.out = i - 1), 1, # control
                                              rep_len(0, length.out = 12 - i), -1, 
                                              rep_len(0, length.out = 12)))
      DGE_glmFTest <- DGE_glmFTest$table
      DGE_glmFTest$FDR <- p.adjust(DGE_glmFTest$PValue,
                                   method = "fdr")
      DGE_glmFTest$Gene_Symbol <- rownames(DGE_glmFTest)
      DGE_glmFTest <- join(DGE_glmFTest,
                           ENSG_coord_gene_gencodev28,
                           by = "Gene_Symbol", 
                           match = "first")
    } else if (i > 13) {
      DGE_glmFTest <- glmQLFTest(DGE_glmFit,
                                 contrast = c(rep_len(0, length.out = 12), -1, # control
                                              rep_len(0, length.out = i - 14), 1, 
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
                file = paste("by_gene/", gRNA_list[i], "_table.txt", sep = ""), 
                row.names = F,
                quote = F, sep = "\t")
  })
}

save.image()
