# make disease/gene group plots
# Siwei 16 Mar 2022

# init
library(ggplot2)
library(readr)
library(RColorBrewer)
library(stringr)

##
raw_data_df <-
  read_delim("df_load.txt", delim = "\t",
             escape_double = FALSE, col_names = FALSE,
             trim_ws = TRUE)
diseases_list <-
  read_csv("disease_order.txt",
           col_names = FALSE)
diseases_list <- unlist(diseases_list)
category_list <-
  unlist(str_split(string = diseases_list,
                   pattern = '_hg19',
                   simplify = T)[ , 1])
df_integrated <- raw_data_df


colnames(df_integrated) <-
  c("Category",
    "Prop._SNPs",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Enrichment_p")

df_integrated$Category <- category_list

deconstructed_names <-
  data.frame(str_split(df_integrated$Category,
                       pattern = "_",
                       simplify = T),
             stringsAsFactors = F)

df_integrated$Disease <- ""
df_integrated$Celltype_timing <- ""

df_integrated$Disease <-
  deconstructed_names$X1
df_integrated$Celltype_timing <-
  apply(deconstructed_names[, 2:5],
        1, paste,
        collapse = "_")
df_integrated$Celltype_timing <-
  str_remove(string = df_integrated$Celltype_timing,
             pattern = "_$")
df_integrated$Celltype_timing <-
  str_remove(string = df_integrated$Celltype_timing,
             pattern = "_$")

# df_integrated$Disease[1:36] <-
#   deconstructed_names$X1[1:36]
# df_integrated$Celltype_timing[1:36] <-
#   apply(deconstructed_names[1:36, 2:6],
#         1, paste,
#         collapse = "_")
#
# df_integrated$Disease[37:48] <-
#   apply(deconstructed_names[37:48, 1:2],
#         1, paste,
#         collapse = "_")
# df_integrated$Celltype_timing[37:48] <-
#   apply(deconstructed_names[37:48, 3:6],
#         1, paste,
#         collapse = "_")
#
# df_integrated$Disease[49:72] <-
#   deconstructed_names$X1[49:72]
# df_integrated$Celltype_timing[49:72] <-
#   apply(deconstructed_names[49:72, 2:4],
#         1, paste,
#         collapse = "_")
#
# df_integrated$Disease[73:84] <-
#   apply(deconstructed_names[73:84, 1:3],
#         1, paste,
#         collapse = "_")
# df_integrated$Celltype_timing[73:84] <-
#   apply(deconstructed_names[73:84, 4:6],
#         1, paste,
#         collapse = "_")
#
# df_integrated$Disease[85:96] <-
#   deconstructed_names$X1[85:96]
# df_integrated$Celltype_timing[85:96] <-
#   apply(deconstructed_names[85:96, 2:4],
#         1, paste,
#         collapse = "_")

###
df_to_plot <- df_integrated
df_to_plot$Enrichment <- abs(df_to_plot$Enrichment)
df_to_plot$`-logP` <- 0 - log10(df_to_plot$Enrichment_p)


df_to_plot$Category <- factor(df_to_plot$Category)
df_to_plot$Disease <- factor(df_to_plot$Disease)
df_to_plot$Celltype_timing <- factor(df_to_plot$Celltype_timing)
# df_to_plot$Contrasts <- unlist(str_split(df_to_plot$VARIABLE,
#                                          pattern = "_sig",
#                                          simplify = T)[, 1])
df_to_plot <-
  df_to_plot[!(df_to_plot$Disease %in% "UC"), ]
df_to_plot$Celltype_timing <-
  str_replace_all(string = df_to_plot$Celltype_timing,
                  pattern = "^NmCp",
                  replacement = "NEFM-")
df_to_plot$Celltype_timing <-
  str_replace_all(string = df_to_plot$Celltype_timing,
                  pattern = "^NpCm",
                  replacement = "NEFM+")

ggplot(df_to_plot,
       aes(x = Celltype_timing,
           y = Disease,
           size = Enrichment,
           fill = `-logP`)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  xlab("Cell_Linktype") +
  scale_y_discrete(limits = rev) +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0)) +
  ggtitle("LDSC Gene-peak link enrichment")

# df_to_plot$ct_test <- df_to_plot$Celltype_timing
# df_to_plot$ct_test <-
#   str_remove(string = df_to_plot$ct_test,
#              pattern = "_$")
