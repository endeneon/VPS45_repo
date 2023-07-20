# Siwei Jan 10 2021
# make ASD BD and SCZ plots to replace Marc's plots
# need to consider whether to use line 11, 12, or merged data

# init
library(ggplot2)
library(RColorBrewer)
library(scales)

# use ASD, BD, and SCZ data from DER_13_Disorder

# make SCZ plot
plot_table <- merge(VPS45_AAvsGG_data,
                    DER_13_Disorder_DEX_Genes_details_DGE,
                    by.x = "Geneid",
                    by.y = "ensembl_gene_id")
plot_table <- plot_table[plot_table$FDR < 0.2, ]
plot_table <- plot_table[plot_table$SCZ.fdr < 0.05, ]

ggplot(plot_table, 
       aes(x = logFC,
           y = SCZ.log2FC)) +
  geom_point(alpha = 0.8,
             size = 0.7,
             colour = brewer.pal(3, "Dark2")[3]) +
  stat_smooth(method = "lm", 
              se = F,
              fullrange = T,
              colour = "black") +
  xlim(-2, 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(-0.4, 0.4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ggtitle("Spearman's rho = 0.313, p < 2.2x10-16") +
  xlab("CROP-seq log2FC") +
  # scale_colour_manual(values = brewer.pal(3, "Dark2")[1]) +
  theme_bw()

cor.test(plot_table$logFC, 
         plot_table$SCZ.log2FC, 
         alternative = "t", 
         method = "s")

# make BD plot
plot_table <- merge(VPS45_AAvsGG_data,
                    DER_13_Disorder_DEX_Genes_details_DGE,
                    by.x = "Geneid",
                    by.y = "ensembl_gene_id")
plot_table <- plot_table[plot_table$FDR < 0.2, ]
plot_table <- plot_table[plot_table$BD.fdr < 0.05, ]

ggplot(plot_table, 
       aes(x = logFC,
           y = BD.log2FC)) +
  geom_point(alpha = 0.8,
             size = 0.7,
             colour = brewer.pal(3, "Dark2")[2]) +
  stat_smooth(method = "lm", 
              se = F,
              fullrange = T,
              colour = "black") +
  xlim(-2, 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(-0.4, 0.4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ggtitle("Spearman's rho = 0.372, p = 5.29x10-7") +
  xlab("CROP-seq log2FC") +
  # scale_colour_manual(values = brewer.pal(3, "Dark2")[1]) +
  theme_bw()

cor.test(plot_table$logFC, 
         plot_table$BD.log2FC, 
         alternative = "t", 
         method = "s")

# make ASD plot
plot_table <- merge(VPS45_AAvsGG_data,
                    DER_13_Disorder_DEX_Genes_details_DGE,
                    by.x = "Geneid",
                    by.y = "ensembl_gene_id")
plot_table <- plot_table[plot_table$FDR < 0.2, ]
plot_table <- plot_table[plot_table$ASD.fdr < 0.05, ]

ggplot(plot_table, 
       aes(x = logFC,
           y = ASD.log2FC)) +
  geom_point(alpha = 0.8,
             size = 0.7,
             colour = brewer.pal(3, "Dark2")[1]) +
  stat_smooth(method = "lm", 
              se = F,
              fullrange = T,
              colour = "black") +
  xlim(-2, 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(-0.4, 0.4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ggtitle("Spearman's rho = 0.377, p = 1.07x10-10") +
  xlab("CROP-seq log2FC") +
  # scale_colour_manual(values = brewer.pal(3, "Dark2")[1]) +
  theme_bw()

cor.test(plot_table$logFC, 
         plot_table$ASD.log2FC, 
         alternative = "t", 
         method = "s")

save.image(file = "plot_FigS7.RData")
