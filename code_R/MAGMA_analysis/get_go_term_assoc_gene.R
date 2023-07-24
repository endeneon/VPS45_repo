# 21 Jun 2019
# Siwei
# get all gene names associated with specific GO term


library(biomaRt)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
View(listAttributes(ensembl))
View(listFilters(ensembl))


#gets gene symbol and go_id for all genes annotated
# filter use go term ('go') and GO term id (values =)
# GO:0007399~nervous system development
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id', 'hgnc_symbol'),
                   filters = 'go', values = 'GO:0007399',
                   mart = ensembl,
                   uniqueRows = T)
gene.data <- gene.data[!duplicated(gene.data$hgnc_symbol), ]
write.table(gene.data[, 3], file = "selected_go_term_associated_genes/GO_0007399_nervous_system_development_with_children_terms.txt",
            quote = F,
            row.names = F, col.names = F,
            sep = "\t")

# GO:0030182~neuron differentiation
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id', 'hgnc_symbol'),
                   filters = 'go', values = 'GO:0030182',
                   mart = ensembl,
                   uniqueRows = T)
gene.data <- gene.data[!duplicated(gene.data$hgnc_symbol), ]
write.table(gene.data[, 3], file = "selected_go_term_associated_genes/GO_0030182_neuron differentiation_with_children_terms.txt",
            quote = F,
            row.names = F, col.names = F,
            sep = "\t")

# GO:0007420~brain development
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id', 'hgnc_symbol'),
                   filters = 'go', values = 'GO:0007420',
                   mart = ensembl,
                   uniqueRows = T)
gene.data <- gene.data[!duplicated(gene.data$hgnc_symbol), ]
write.table(gene.data[, 3], file = "selected_go_term_associated_genes/GO_0007420_brain development_with_children_terms.txt",
            quote = F,
            row.names = F, col.names = F,
            sep = "\t")

