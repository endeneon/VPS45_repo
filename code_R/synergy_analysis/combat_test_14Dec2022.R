# 14 Dec 2022 Siwei
# test ComBat function with reference batch

library(sva)
library(bladderbatch)

data(bladderdata)

dat <- bladderEset[1:50,]

pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer),
                   data = pheno)
# reference-batch version, with covariates
combat_edata3 =
  ComBat(dat = edata,
         batch = batch,
         mod = mod,
         par.prior = F,
         ref.batch = 3)
