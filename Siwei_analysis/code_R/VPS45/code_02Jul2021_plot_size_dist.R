# Siwei 02 Jul 2021
# plot read size distribution of the 20 lines

# init
library(readr)
library(ggplot2)
library(RColorBrewer)

library(stringr)


# load data
size.data.files <- 
  list.files(path = "size_dist/",
             pattern = "*.txt")

list.size.data <- vector(mode = "list", length = length(size.data.files))

i <- 1
for (i in 1:length(size.data.files)) {
  list.size.data[[i]] <-
    read_delim(paste0("size_dist/",
                      size.data.files[i]),  
               delim = "\t", 
               escape_double = FALSE, 
               trim_ws = TRUE, 
               skip = 10)
}
names(list.size.data) <- 
  str_replace_all(str_remove_all(size.data.files, pattern = "\\.txt"),
                 pattern = "Glut_rapid_neuron20_", 
                 replacement = "Cell line ")

# add percentage to list.size.data
i <- 1
for (i in 1:length(size.data.files)) {
  list.size.data[[i]]$percentage <-
    list.size.data[[i]]$All_Reads.fr_count * 100 / sum(list.size.data[[i]]$All_Reads.fr_count)
  
  list.size.data[[i]]$max.percent <-
    list.size.data[[i]]$All_Reads.fr_count * 100 / max(list.size.data[[i]]$All_Reads.fr_count)

}

# make plot
colour.list <- RColorBrewer::brewer.pal(10, "Paired")

##
i <- 1
for (i in 1:10) {
  if (i == 1) {
    plot.size.dist <- 
      ggplot() +
      geom_line(data = list.size.data[[i]],
                aes(x = insert_size,
                    y = max.percent),
                colour = colour.list[i],
                alpha = 0.7)
  } else {
    plot.size.dist <- 
      plot.size.dist +
      # ggplot() + # note no new ggplot() is required here
      geom_line(data = list.size.data[[i]],
                aes(x = insert_size,
                    y = max.percent),
                colour = colour.list[i],
                alpha = 0.7)

  }
}

plot.size.dist <-
  plot.size.dist +
  theme_classic() +
  xlab("Insert size") + 
  ylab(("Percentage")) +
  xlim(0, 800)
  
plot.size.dist

##
rm(plot.size.dist)
i <- 11
for (i in 11:20) {
  print(i)
  if (i == 11) {
    plot.size.dist <- 
      ggplot() +
      geom_line(data = list.size.data[[i]],
                aes(x = insert_size,
                    y = max.percent),
                colour = colour.list[i - 10],
                alpha = 0.7)
  } else {
    plot.size.dist <- 
      plot.size.dist +
      # ggplot() + # note no new ggplot() is required here
      geom_line(data = list.size.data[[i]],
                aes(x = insert_size,
                    y = max.percent),
                colour = colour.list[i - 10],
                alpha = 0.7)
    
  }
}

plot.size.dist <-
  plot.size.dist +
  theme_classic() +
  xlab("Percentage") + 
  ylab(("Read count")) +
  xlim(0, 800)

plot.size.dist
