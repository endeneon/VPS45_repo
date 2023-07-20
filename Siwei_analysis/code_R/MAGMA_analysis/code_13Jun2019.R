# 13 Jun 2019
# test magma.celltyping

library(MAGMA.Celltyping)
data("ctd")
data("ortholog_data_Mouse_Human")

print(ctd[[1]]$quantiles[c("Gfap","Dlg4","Aif1"),])


View(PGC_SCZ_full_hg19[1:10,])

write.table(PGC_SCZ_full_hg19[, c(1, 13)],
            file = "PGC_SCZ_hg19_pval.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")
