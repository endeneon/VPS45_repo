# plot peak size distribution of iN-Glut 20
# Siwei 09 Jul 2021

# init
library(readr)
library(stringr)

# load data
df.raw <- read_delim("size_dist/NGN2_Glut_20_lines_21Jun2021_peaks.narrowPeak", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
df.raw <- 
  df.raw[str_detect(string = df.raw$X1, # remove peaks do not start with chr
                    pattern = "^chr[:digit:]"), ]
df.raw <-
  df.raw[!(str_detect(string = df.raw$X1, # remove peaks in alt assembly
                    pattern = "_")), ]
unique(df.raw$X1)

df.raw$peak.width <- df.raw$X3 - df.raw$X2

hist(df.raw$peak.width, 
     breaks = 1000000, 
     include.lowest = T, 
     xlim = c(0, 3000),
     xlab = "Peak width (bp)", 
     ylab = "Counts")

mean(df.raw$peak.width)
median(df.raw$peak.width)


## calculate CN peaks
# load data
df.raw <- read_delim("size_dist/CN_all_peaks.narrowPeak", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
df.raw <- 
  df.raw[str_detect(string = df.raw$X1, # remove peaks do not start with chr
                    pattern = "^chr[:digit:]"), ]
df.raw <-
  df.raw[!(str_detect(string = df.raw$X1, # remove peaks in alt assembly
                      pattern = "_")), ]
unique(df.raw$X1)

df.raw$peak.width <- df.raw$X3 - df.raw$X2

hist(df.raw$peak.width, 
     breaks = 1000000, 
     include.lowest = T, 
     xlim = c(0, 3000),
     xlab = "Peak width (bp)", 
     ylab = "Counts")

mean(df.raw$peak.width)
median(df.raw$peak.width)
