
library(tidyverse)
library(ggplot2)
genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_csv2 <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv")
genedata_df <- t(data.frame(genedata_csv))
genedata_df2 <- t(data.frame(genedata_csv2))
genedata_matrix <- data.matrix(genedata_df)
genedata_means <- read.csv("Means.csv")
genedata_stdevs <- read.csv("Standard Deviations.csv")
genedata_genes <- genedata_df2[c(2:19178),1]
dupes <- c()

for (gene in genedata_genes){
  if (is.null(dupes[gene])){
    dupes[gene] <- 1
  }
  else {
    dupes[gene] <- dupes[gene] + 1
  }
}

dupes

if(is.null(dupes[genedata_genes[2]])){
  dupes[genedata_genes[2]] <- 1
}
dupes

genedata_genes[2]

