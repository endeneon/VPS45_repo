

test_table <- merge(VPS45_all_lines, 
                    RNAseq_ASD_4region_sumstats, 
                    by.x = "Geneid", by.y = "X1")
test_table$logFC.groupGG <- 0 - test_table$logFC.groupGG
test_table <- test_table[test_table$Frontal.P.Value < 0.05, ]
test_table <- test_table[test_table$FDR < 0.05, ]

cor.test(test_table$logFC.groupGG, 
         test_table$Frontal.logFC, 
         alternative = "t", 
         method = "s")
0 - log10(cor.test(test_table$logFC.groupGG, 
                   test_table$Frontal.logFC, 
                   alternative = "t", 
                   method = "s"))


ggplot(test_table, 
        aes(x = logFC,
            y = Frontal.logFC)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

test_table <- merge(VPS45_all_lines, 
                    DER_13_Disorder_DEX_Genes_details_DGE, 
                    by.x = "Geneid", 
                    by.y = "ensembl_gene_id", no.dups = T)

# test_table <- test_table[test_table$Frontal.P.Value < 0.05, ]

test_table$FDR_SCZ <- p.adjust(test_table$SCZ.p.value,
                               method = "fdr")
test_table$logFC.groupGG <- 0 - test_table$logFC.groupGG
test_table <- test_table[test_table$FDR < 0.05, ]
test_table <- test_table[test_table$FDR_SCZ < 0.05, ]

cor.test(test_table$logFC.groupGG, 
         test_table$SCZ.log2FC, 
         alternative = "t", 
         method = "s")

0 - log10(cor.test(test_table$logFC.groupGG, 
                   test_table$SCZ.log2FC, 
                   alternative = "t", 
                   method = "s")$p.value)

ggplot(test_table, 
       aes(x = logFC.groupGG,
           y = SCZ.log2FC)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

plot_table_SCZ <- test_table[, c("logFC.groupGG", "SCZ.log2FC")]
plot_table_SCZ$Disease <- "SCZ"

write.table(test_table, 
            file = "SCZ_FDR005_382_genes.txt", 
            row.names = F, col.names = T, 
            sep = "\t", quote = F)

###

test_table <- merge(VPS45_all_lines, 
                    DER_13_Disorder_DEX_Genes_details_DGE, 
                    by.x = "Geneid", 
                    by.y = "ensembl_gene_id", no.dups = T)
test_table$FDR_BD <- p.adjust(test_table$BD.p.value,
                               method = "fdr")
test_table$logFC.groupGG <- 0 - test_table$logFC.groupGG
test_table <- test_table[test_table$FDR < 0.05, ]
test_table <- test_table[test_table$FDR_BD < 0.05, ]

cor.test(test_table$logFC.groupGG, 
         test_table$BD.log2FC, 
         alternative = "t", 
         method = "s")

0 - log10(cor.test(test_table$logFC.groupGG, 
                   test_table$BD.log2FC, 
                   alternative = "t", 
                   method = "s")$p.value)

plot_table_BP <- test_table[, c("logFC.groupGG", "BD.log2FC")]
plot_table_BP$Disease <- "Bip"

ggplot(test_table, 
       aes(x = logFC.groupGG,
           y = BD.log2FC)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

write.table(test_table, 
            file = "BP_FDR005_113_genes.txt", 
            row.names = F, col.names = T, 
            sep = "\t", quote = F)

###

test_table <- merge(VPS45_all_lines, 
                    DER_13_Disorder_DEX_Genes_details_DGE, 
                    by.x = "Geneid", 
                    by.y = "ensembl_gene_id", no.dups = T)
test_table$FDR_ASD <- p.adjust(test_table$ASD.p.value,
                              method = "fdr")
test_table$logFC.groupGG <- 0 - test_table$logFC.groupGG
test_table <- test_table[test_table$FDR < 0.05, ]
test_table <- test_table[test_table$FDR_ASD < 0.05, ]

cor.test(test_table$logFC.groupGG, 
         test_table$ASD.log2FC, 
         alternative = "t", 
         method = "s")

0 - log10(cor.test(test_table$logFC.groupGG, 
                   test_table$ASD.log2FC, 
                   alternative = "t", 
                   method = "s")$p.value)

plot_table_ASD <- test_table[, c("logFC.groupGG", "ASD.log2FC")]
plot_table_ASD$Disease <- "ASD"

ggplot(test_table, 
       aes(x = logFC.groupGG,
           y = ASD.log2FC)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

write.table(test_table, 
            file = "ASD_FDR005_201_genes.txt", 
            row.names = F, col.names = T, 
            sep = "\t", quote = F)

###


save.image(file = "scatter_plot.RData")


colnames(plot_table_ASD) <- c("logFC.VPS45", "logFC.Disease", "Disease")
colnames(plot_table_BP) <- c("logFC.VPS45", "logFC.Disease", "Disease")
colnames(plot_table_SCZ) <- c("logFC.VPS45", "logFC.Disease", "Disease")

plot_table <- data.frame(rbind(plot_table_ASD,
                               plot_table_BP,
                               plot_table_SCZ),
                         stringsAsFactors = F)

ggplot(plot_table, 
       aes(x = logFC.VPS45,
           y = logFC.Disease,
           colour = Disease)) +
  geom_point(alpha = 0.8,
             size = 0.7) +
  stat_smooth(method = "lm", 
              se = F,
              fullrange = T) +
  xlim(-2, 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(-0.75, 0.75) +
  # geom_abline(slope = 1, 
  #             intercept = 0, linetype = "dashed") +
  scale_colour_manual(labels = c("ASD, rho=0.298, p=1.91e-5, n=201",
                                 "Bip, rho=0.371, p=5.74e-5, n=113",
                                 "SCZ, rho=0.336, p=2.04e-11, n=382"),
                      values = brewer.pal(3, "Dark2")) +
  # scale_colour_discrete(labels = c("ASD, rho=0.298, p=1.91e-5, n=201",
  #                                  "Bip, rho=0.371, p=5.74e-5, n=113",
  #                                  "SCZ, rho=0.336, p=2.04e-11, n=382"),
  #                       values = brewer.pal(3, "Dark2")) +
  theme_bw() +
  theme(legend.position = c(0.25, 0.9), 
        legend.direction = "vertical",
        axis.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha(colour = "white",
                                                      alpha = 0.7)))

ggplot(plot_table, 
       aes(x = logFC.VPS45,
           y = logFC.Disease,
           colour = Disease)) +
  geom_point(alpha = 0.8,
             size = 0.7) +
  stat_smooth(method = "lm", 
              se = F,
              fullrange = T) +
  xlim(-2, 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(-0.75, 0.75) +
  # geom_abline(slope = 1, 
  #             intercept = 0, linetype = "dashed") +
  scale_colour_manual(labels = c("ASD, rho=0.30, p=1.9e-5",
                                 "BD, rho=0.37, p=5.7e-5",
                                 "SZ, rho=0.34, p=2.0e-11"),
                      values = brewer.pal(3, "Dark2")) +
  # scale_colour_discrete(labels = c("ASD, rho=0.298, p=1.91e-5, n=201",
  #                                  "Bip, rho=0.371, p=5.74e-5, n=113",
  #                                  "SCZ, rho=0.336, p=2.04e-11, n=382"),
  #                       values = brewer.pal(3, "Dark2")) +
  theme_bw() +
  theme(legend.position = c(0.21, 0.88), 
        legend.direction = "vertical",
        axis.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha(colour = "white",
                                                      alpha = 0.7)))


RColorBrewer::brewer.pal.info
brewer.pal(3, "Dark2")
