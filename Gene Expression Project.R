genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_df <- (data.frame(genedata_csv))

genedata_df2 <- t(data.frame(genedata_csv))
genedata_genes <- genedata_df[1,]
head(genedata_genes, 20)

library(tidyverse)

genedata_df[nrow(genedata_df)+1,] <- genedata_stdevs
genedata_df[nrow(genedata_df),1393] <- "stdev"

genedata_stdevs <- apply(genedata_df[,-1], 1, sd, na.rm = TRUE)



genedata_stdevs

genedata_stdevs_df <- data.frame(genedata_stdevs)


library(ggplot2)

install.packages("tidyverse")
install.packages("ggplot2")

ggplot(data=NULL, aes(x=genedata_stdevs)) + geom_histogram(bins = 200)
