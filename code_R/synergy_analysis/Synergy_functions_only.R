
## 3.	Load R packages

pacman::p_load(limma, edgeR, pheatmap, RColorBrewer,
               ggplot2, ggpubr, qvalue,  plyr, wesanderson,
               GSEABase, grid, scales, WebGestaltR, stringr)


## 4.	Run custom functions

# The mds() function is based on plotMDS() in the limma package.
# When given a DGE object and a column in the meta data table, containing groups of interest,
# it produces an multidimensional scaling plot, colored by the provided groups.

mds <- function(normDGE, metacol, title){
  mcol <- as.factor(metacol)
  col <- rainbow(length(levels(mcol)), 1, 0.8, alpha = 0.5)[mcol]
  plotMDS(normDGE, col = col, pch = 16, cex = 2)
  legend("center",
         fill = rainbow(length(levels(mcol)), 1, 0.8),
         legend = levels(mcol),
         horiz = F,
         bty = "o",
         box.col="grey",
         xpd=TRUE)
  title(main=title)
}


# The cameraplusplots() function is based on camera() in the limma package.
# When given a contrast in form of a named vector (such as a column of the contrast matrix),
#a list of gene set groups, a voom object, a design matrix and a color palette for the gene set groups,
# it produces a scatter plot of all tested gene sets and their adjusted P-values
# as well as a bar graph of the 10 most significant gene sets, colored by gene set group/category.

cameraplusplots <- function(contrast, genesetlist, vobject, design, catcolors, title){
  tmp.list <- list()
  cam <- data.frame(matrix(ncol = 5, nrow = 0))
  for (i in 1:length(genesetlist)){
    cam.s <- camera(vobject, genesetlist[[i]], design, contrast = contrast, inter.gene.cor = 0.01)
    tmp.list[[i]] <- cam.s
    names(tmp.list)[i] <- names(genesetlist)[i]
    tmp.list[[i]]$category <- names(tmp.list[i])
    colnames(cam) <- names(tmp.list[[1]])
    cam <- rbind.data.frame(cam, tmp.list[[i]])
    print(paste0("Gene set categories run: ", i))
  }
  cam$neglogFDR <- -log10(cam$FDR)
  ## for plotting purposes only:
  cam$dirNeglogFDR <- cam$neglogFDR
  cam[(cam$Direction == "Down"), "dirNeglogFDR"] <-
    -cam[(cam$Direction == "Down"), "neglogFDR"]
  grob <-
    grobTree(textGrob(c("UP","DOWN"),
                      x = c(0.94, 0.89),
                      y = c(0.95, 0.05),
                      hjust = 0,
                      gp = gpar(fontsize = 13)))
  q <-
    ggplot(aes(x = cam$category,
               y = dirNeglogFDR,
               color = category),
           data = cam) +
    scale_color_manual(values = catcolors) +
    geom_jitter(aes(size = NGenes,
                    alpha = neglogFDR),
                pch = 19,
                show.legend = F) +
    scale_size_continuous(range = c(4,16)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    geom_hline(yintercept = c(-1.3, 1.3),
               color = "red",
               alpha = 0.5) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-10, 10),
                       oob = squish,
                       labels = abs) +
    labs(x = "Gene set categories", y = "-log10(FDR)", title = title) +
    theme_bw(14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    annotation_custom(grob)
  print(q)
  cam$geneSet <- row.names(cam)
  cam10 <- as.data.frame(cam)
  cam10 <- cam10[order(cam10$FDR),]
  cam10 <- cam10[1:10,]
  grob <-
    grobTree(textGrob(c("DOWN","UP"),
                      x = c(0.03, 0.9),
                      y=c(0.025), hjust = 0,
                      gp = gpar(fontsize = 9,
                                col = "grey60")))
  g <-
    ggplot(aes(x = geneSet,
               y = dirNeglogFDR,
               fill = category),
           data = cam10) +
    geom_col() +
    aes(reorder(stringr::str_wrap(geneSet, 60),-FDR), dirNeglogFDR) +
    xlab(NULL) +
    geom_hline(yintercept = c(-1.3, 1.3), color = "red", alpha = 0.3) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-10, 10), oob = squish, labels = abs) +
    labs(y = "-log10(FDR)", title = title) +
    scale_fill_manual(values = catcolors) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    annotation_custom(grob)
  print(g)
  return(cam)
}


# The oraplot() function takes the data frame resulting from over-representation analysis using WebGestaltR
# and a color palette as input and returns a bar graph of the 10 most significant gene sets.

oraplot <- function(orares, catcolors, name){
  orares.n <- orares[order(orares$FDR), ]
  orares.n <- orares.n[1:10, ]
  orares.n$neglogFDR <- -log10(orares.n$FDR)
  orares.n <- orares.n[orares.n$neglogFDR>0, ]
  orares.n$geneSet <- gsub("_", " ", orares.n$geneSet)
  g <-
    ggplot(aes(x=reorder(str_wrap(geneSet, 60), neglogFDR),
               y = neglogFDR,
               fill = category),
           data = orares.n) +
    geom_col() +
    geom_hline(yintercept = 1.3, color = "red", alpha = 0.5) +
    labs(y = "-log10(FDR)", x = "", title = paste0(name)) +
    scale_fill_manual(values = catcolors) +
    coord_flip() +
    theme_bw(11)
  return(g)
}


# The power.compare.logFC() function is given the variances of the combinatorial perturbation
# and the additive model comparisons (sig1, sig2; variance in the additive model is usually higher,
# in proportion to the number of individual perturbations).
# Further, the number of samples used to determine the variances (N), a vector of sample numbers of interest (N_other),
# a significance cutoff (alpha) and the number of tests performed (n_tests, usually the number of transcripts).
# N_other/N is the relative sample size and the variance of logFC1 - logFC2 is sig1^2 + sig2^2.
# Since the standard error of the mean is inversely proportional to sqrt(N), multiplying the sample size by F decreases the SE by sqrt(F).
# On the variance scale, this corresponds to dividing by n_scale.

power.compare.logFC <- function(sig1, sig2, N, N_other = c(2,4,6,8,10), alpha = 0.05, n_tests = 20000){
  d <- seq(0, 3, length.out=1000)
  alpha_multiple <- alpha / n_tests
  df <- lapply( N_other/N, function(n_scale){
    sigSq <- (sig1^2 + sig2^2) / n_scale
    cutoff <- qnorm( alpha_multiple/2, 0, sd = sqrt(sigSq), lower.tail = FALSE)
    p1 <- pnorm(-1*cutoff, d, sqrt(sigSq))
    p2 <- 1-pnorm(cutoff, d, sqrt(sigSq))
    data.frame(n_scale, d, power=p1+p2)
  })
  df <- do.call("rbind", df)
  ggplot(df, aes(d, power, color = as.factor(n_scale*N))) +
    geom_line() +
    theme_bw(14) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ylim(0, 1) +
    scale_color_discrete("Samples") +
    xlab(bquote(abs(logFC[observed] - logFC[expected]))) +
    ggtitle("Power versus difference in logFC")
}


# The categorize.synergy() function is given the combined matrix of log2FC values of the additive model,
# the combinatorial perturbation and the synergistic effect differential expression results.
# It creates a new column in the resulting data frame, assigning synergy categories to each gene.

categorize.synergy <- function(logFCmatrix, meanSE){
  m <- logFCmatrix
  m$magnitude.syn <- NA
  for (i in 1:length(m$Gene_name)){
    if (m$Synergistic.logFC[i] > meanSE){
      if (m$Additive.logFC[i] < -meanSE){
        if (m$Combinatorial.logFC[i] > meanSE){
          m$magnitude.syn[i] = "more.up"
        } else m$magnitude.syn[i] = "less.down"
      } else m$magnitude.syn[i] = "more.up"
    }
    else if (m$Synergistic.logFC[i] < -meanSE){
      if (m$Additive.logFC[i] > meanSE){
        if (m$Combinatorial.logFC[i] < -meanSE){
          m$magnitude.syn[i] = "more.down"
        } else m$magnitude.syn[i] = "less.up"
      } else m$magnitude.syn[i] = "more.down"
    } else m$magnitude.syn[i] = "same"
  }
  m$magnitude.syn <- as.factor(m$magnitude.syn)
  return(m)
}


# The stratify.by.syn.cat() function is given a subset of interest of the table created by the categorize.synergy() function.
# It creates a list object containing vectors of gene symbols by synergy category,
# which is used as input for over-representation analysis with WebGestaltR.

stratify.by.syn.cat <- function(log2FC.matrix.sub){
  synergy.cat.list <- list("less.down" = as.character(log2FC.matrix.sub[
    log2FC.matrix.sub$magnitude.syn == "less.down", "Gene_name"]),
    "less.up" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "less.up", "Gene_name"]),
    "more.down" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "more.down", "Gene_name"]),
    "more.up" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "more.up", "Gene_name"]),
    "same" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "same", "Gene_name"]))
  return(synergy.cat.list)
}


