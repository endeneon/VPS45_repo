# 21 Mar 2022
# Siwei make ldsc gene plots

# init
library(readr)
library(ggplot2)
library(RColorBrewer)
library(stringr)

result_files <- dir(path = "ldsc_results/",
                    pattern = "results$")

raw_df_list <- vector(mode = "list",
                      length = length(result_files))
names(raw_df_list) <- result_files
# read in all data
i <- 1
for (i in 1:length(result_files)) {
  raw_df_list[[i]] <-
    read_delim(paste0('ldsc_results/',result_files[i]), 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE)
}

# extract row 2
i <- 1
for (i in 1:length(result_files)) {
  if (i == 1) {
    raw_df_aggregated <-
      raw_df_list[[i]][2, ]
  } else {
    raw_df_aggregated <-
      rbind(raw_df_aggregated,
            raw_df_list[[i]][2, ])
  }
}
# rownames(raw_df_aggregated) <- result_files
raw_df_aggregated$Category <- result_files
raw_df_aggregated$'-log10P' <-
  0 - log10(raw_df_aggregated$Enrichment_p)

raw_df_aggregated$Disease <- 
  unlist(str_split(string = raw_df_aggregated$Category,
                   pattern = "_",
                   simplify = T)[, 1])
# raw_df_aggregated$GeneSet <- 
#   paste(str_split(string = raw_df_aggregated$Category,
#                    pattern = "_",
#                    n = Inf,
#                    simplify = T)[, 2:4],
#         collapse = "_")
raw_df_aggregated$GeneSet <-
  rep_len(str_sub(string = raw_df_aggregated$Category[1:16],
                  start = 6L),
          length.out = 112)

raw_df_aggregated$Disease[raw_df_aggregated$Disease %in% "daner"] <- "daner_SCZ3"
raw_df_aggregated$Disease[raw_df_aggregated$Disease %in% "PGC"] <- "PGC_ASD"

ggplot(raw_df_aggregated,
       aes(x = Disease,
           y = GeneSet,
           size = Enrichment,
           fill = `-log10P`)) +
  geom_point(shape = 21) +
  scale_fill_gradient2(low = "white",
                       high = "red") +
  theme_bw()
