# Siwei 09 Dec 2022
# verify C1orf54 transcript expression by Kallisto


# init
library(readr)
library(readxl)
library(tximport)

library(stringr)
library(edgeR)
library(variancePartition)
library(factoextra)

library(ggplot2)

library(GenomicFeatures)

# get the list of all kallisto output files
input_file_list <-
  list.files(path = "/data/FASTQ/RNASeq_Novogene_05Dec2022/Novogene_VPS45_3xKD_05Dec2022/01.RawData/kallisto_output/",
             pattern = "abundance.tsv",
             full.names = T,
             recursive = T)

names(input_file_list) <-
  str_split(string = input_file_list,
            pattern = "\\/",
            simplify = T)[, 9]

df_from_kallisto <-
  tximport(files = input_file_list,
           type = "kallisto",
           txOut = T,
           ignoreAfterBar = T)

df_raw <-
  df_from_kallisto$counts
df_raw <-
  as.data.frame(df_raw)

df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST00000369099"), ] # 202, 509 bp
df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST00000369102.*"), ] # 203, 1251 bp
df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST00000369098.*"), ] # 201, 418 bp

df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST0000036913.*"), ] # 
df_raw <-
  df_from_kallisto$abundance
df_raw <-
  as.data.frame(df_raw)

df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST00000369099"), ] # 202, 509 bp
df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST00000369102.*"), ] # 203, 1251 bp
df_raw[str_detect(string = rownames(df_raw),
                  pattern = "ENST00000369098.*"), ] # 201, 418 bp
