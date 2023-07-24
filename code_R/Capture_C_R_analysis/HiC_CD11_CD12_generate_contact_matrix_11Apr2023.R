# Siwei 11 Apr 2023
# analyse Capture C and CHiCAGO output of CD11/CD12 output HiC data

library(Chicago)
library(plotgardener)

library("stringr")

# library(remotes)
# remotes::install_github("PhanstielLab/BentoBox")

strawr::readHicBpResolutions("../Capture/hicFiles/CD11-0.hic")

rs2027349_region <-
  readHic(file = "../Capture/hicFiles/CD11-0.hic",
          chrom = "chr1",
          chromstart = 149000000,
          chromend = 151000000,
          assembly = "hg38",
          resolution = 5000,
          matrix = "log2oe",
          norm = "KR")

sub_count_region <-
  rs2027349_region[!(rs2027349_region$chr1_A == rs2027349_region$chr1_B), ]

# hist(sub_count_region$counts)
# sum(sub_count_region$counts > 3)

sub_count_region <-
  sub_count_region[sub_count_region$counts > 3, ]
sub_count_region <-
  sub_count_region[abs(sub_count_region$chr1_B - sub_count_region$chr1_A) > 11000, ]



writeout_df <-
  sub_count_region
writeout_df$chr <- "chr1"
writeout_df <-
  writeout_df[, c(4, 1, 2, 3)]

writeout_df$chr1_A1 <-
  writeout_df$chr1_A + 1

writeout_df$chr1_B <-
  str_c("chr1",
        ":",
        writeout_df$chr1_B,
        "-",
        writeout_df$chr1_B + 1,
        ",",
        writeout_df$counts)

writeout_df <-
  writeout_df[, c(1, 2, 5, 3, 4)]
writeout_df$chr1_A1 <-
  as.character(writeout_df$chr1_A1)

writeout_df$strand <- "."
# writeout_df$chr1_A <-
#   str_c("chr1",
#         writeout_df$chr1_A,
#         writeout_df$chr1_A + 1,
#         sep = ",")
# writeout_df$chr1_B <-
#   str_c("chr1",
#         writeout_df$chr1_B,
#         writeout_df$chr1_B + 1,
#         sep = ",")

write.table(writeout_df,
            file = "rs2027349_long_range_new_browser.txt",
            row.names = F, col.names = F,
            quote = F, sep = "\t")
