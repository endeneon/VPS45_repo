# Siwei 30 Sept 2020
# make the bar plot

# init
library(ggplot2)
library(scales)

##
colnames(table_for_correlation_plotting)[4] <- 'genes'
colnames(table_for_correlation_plotting)[3] <- '-log10P'
colnames(table_for_correlation_plotting)[1] <- 'Diseases'
colnames(table_for_correlation_plotting)[6] <- '-log10P (FDR<0.05)'
colnames(table_for_correlation_plotting)[9] <- "pVal (FDR<0.05)"

table_for_correlation_plotting$pVal <- 
  10^(0 - table_for_correlation_plotting$`-log10P`)
table_for_correlation_plotting$`pVal (FDR<0.05` <- 
  10^(0 - table_for_correlation_plotting$`-log10P (FDR<0.05)`)

table_for_correlation_plotting$Diseases <-
  factor(table_for_correlation_plotting$Diseases,
         levels = c("SCZ_PsychENCODE+GTEx", "SCZ_Fetal_Brain", "Autism_PsychENCODE+GTEx",
                    "Bipolar_PsychENCODE+GTEx", "Alcohol_Dependence", "Major_Depression", 
                    "IBD"))

## make plot for non-FDR trimming
ggplot(table_for_correlation_plotting,
       aes(x = Diseases,
           y = `Spemann's Rho`)) +
  geom_bar(stat = "identity",
           fill = "darkblue", 
           position = position_dodge()) +
  geom_text(aes(label = paste("P=", format(pVal, 
                                           digits = 2,
                                           trim = T),
                              ", N=", format(genes,
                                           digits = 2,
                                           trim = T),
                              sep = "")),
            colour = "black",
            position = position_dodge(1),
            hjust = -0.1,
            size = 3.5) +
  theme_classic() +
  ylim(c(0, 0.25)) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev(levels(table_for_correlation_plotting$Diseases))) +
  coord_flip()

## make plot with FDR trimming
ggplot(table_for_correlation_plotting,
       aes(x = Diseases,
           y = `Spemann's Rho (FDR<0.05)`)) +
  geom_bar(stat = "identity",
           fill = "darkblue", 
           position = position_dodge()) +
  geom_text(aes(label = paste("P=", format(`pVal (FDR<0.05)`, 
                                           digits = 2,
                                           trim = T),
                              ", N=", format(`genes (FDR<0.05)`,
                                             digits = 2,
                                             trim = T),
                              sep = "")),
            colour = "black",
            position = position_dodge(1),
            hjust = -0.1,
            size = 3.5) +
  theme_classic() +
  ylim(c(0, 0.6)) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev(levels(table_for_correlation_plotting$Diseases))) +
  coord_flip()
