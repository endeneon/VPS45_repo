# Siwei 16 Dec 2022
# Run correlation analysis using data in Gandal MJ Science 2018
# and method in Schrode Nat Genet 2019
# calc correlation

## analyse EtOH Alc dependence GSE29555
library(illuminaio)
library(GEOquery)

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE29555", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11001100101011001100101010100100XXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Alc","ctrl"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

# calculate precision weights and show plot of mean-variance trend
# v <- vooma(gset, design, plot=T)
v <- voom(gset, design, plot=T)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.1)

# Collect the information of all genes for correlation analysis
tT <-
  topTable(fit2,
           adjust.method = "fdr",
           sort.by = "B",
           number = length(fit2$F))
tT <-
  subset(tT,
         select = c("ID",
                    "ILMN_Gene",
                    "Symbol",
                    "adj.P.Val",
                    "P.Value",
                    "logFC"))

tT <-
  tT[tT$Symbol != "", ]
sum(tT$adj.P.Val < 0.05)

save(list = c("tT"),
     file = "df_logFC_Alc_Dependence.RData")

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GI","SEQUENCE","GB_ACC"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE29555", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE29555", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 13, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=13", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


### COS data lnalysis
library(readr)
library(stringr)

library(edgeR)

df_COS_raw <-
  read_csv("TF_libraries/GSE106589_geneCounts.csv.gz")
genelist_df_COS_raw <-
  df_COS_raw$gene
df_COS_raw$gene <- NULL
rownames(df_COS_raw) <-
  genelist_df_COS_raw

GSE106589_series_matrix <-
  read_delim("TF_libraries/GSE106589_series_matrix.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE, skip = 27)

df_metadata_COS <-
  as.data.frame(t(GSE106589_series_matrix))
title_df_metadata <- GSE106589_series_matrix$`!Sample_title`
title_df_metadata <-
  str_remove(string = title_df_metadata,
             pattern = "^\\!")
df_metadata_COS <-
  df_metadata_COS[-1, ]
colnames(df_metadata_COS) <-
  title_df_metadata
df_metadata_COS$Sample_name <-
  rownames(df_metadata_COS)

# df_metadata_COS:
# col 9: ctrl/COS
# col 10: sex
# col 11: sample type (NPC/fetal brain)
colnames(df_metadata_COS)[9] <- "COS_phenotype"
colnames(df_metadata_COS)[10] <- "indiv_sex"
colnames(df_metadata_COS)[11] <- "Sample_type_fetal_NPC"

# remove all NPC samples
df_metadata_COS <-
  df_metadata_COS[(df_metadata_COS$Sample_type_fetal_NPC %in% "cell type: 6 wk FB neuron"), ]

df_COS_raw <-
  df_COS_raw[, colnames(df_COS_raw) %in% df_metadata_COS$Sample_name]

# make the raw data and metadata col consistent
df_metadata_COS <-
  df_metadata_COS[df_metadata_COS$Sample_name %in% colnames(df_COS_raw), ]
df_COS_raw <-
  df_COS_raw[,
             match(df_metadata_COS$Sample_name,
                   colnames(df_COS_raw))]
rownames(df_COS_raw) <-
  genelist_df_COS_raw

# replace geneid with gene names
gencode_v35_ENSG_Genename_final <-
  read_delim("gencode.v35.ENSG.Genename.final.list",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[, c(1, 5)]
colnames(gencode_v35_ENSG_Genename_final) <-
  c("Geneid", "Gene_symbol")


gencode_v35_ENSG_Genename_final$Geneid <-
  str_split(gencode_v35_ENSG_Genename_final$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Geneid), ]
gencode_v35_ENSG_Genename_final <-
  gencode_v35_ENSG_Genename_final[!duplicated(gencode_v35_ENSG_Genename_final$Gene_symbol), ]

df_gencode_COS_table <-
  gencode_v35_ENSG_Genename_final[gencode_v35_ENSG_Genename_final$Geneid %in%
                                    rownames(df_COS_raw), ]

df_COS_4_DGE <-
  df_COS_raw[rownames(df_COS_raw) %in% df_gencode_COS_table$Geneid, ]
rownames(df_COS_4_DGE) <-
  genelist_df_COS_raw[rownames(df_COS_raw) %in% df_gencode_COS_table$Geneid]

df_gencode_COS_table <-
  df_gencode_COS_table[match(rownames(df_COS_4_DGE),
                             df_gencode_COS_table$Geneid),
                       ]

rownames(df_COS_4_DGE) <-
  df_gencode_COS_table$Gene_symbol



df_DGE <-
  DGEList(counts = as.matrix(df_COS_4_DGE),
          samples = colnames(df_COS_4_DGE),
          group = unlist(df_metadata_COS[, 9]),
          genes = rownames(df_COS_4_DGE),
          remove.zeros = T)


df_design <-
  model.matrix(~ 0 +
                 factor(indiv_sex) +
                 factor(COS_phenotype),
               data = df_metadata_COS)

df_DGE <-
  calcNormFactors(df_DGE)
df_DGE <-
  estimateDisp(df_DGE,
               design = df_design,
               robust = T)

result_QLM <-
  glmQLFit(df_DGE, design = df_design)
result_QLM <-
  glmQLFTest(result_QLM,
             coef = 3)
results_table <-
  result_QLM$table
results_table$FDR <-
  p.adjust(results_table$PValue,
           method = "fdr")
results_table_COS <- results_table
save(list = c("results_table_COS"),
     file = "df_results_table_COS.RData")


#### CMC data analysis
####


