
library(tidyverse)
library(ggplot2)
genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_df <- t(data.frame(genedata_csv))
genedata_matrix <- data.matrix(genedata_df)
genedata_means <- read.csv("Means.csv")
genedata_stdevs <- read.csv("Standard Deviations.csv")
genedata_genes <- genedata_df[1,]

