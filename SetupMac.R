
library(tidyverse)
library(ggplot2)
genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_csv2 <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv")
genedata_means <- read.csv("Means.csv")[2]
genedata_stdevs <- read.csv("Standard Deviations.csv")[2]
genedata_genes <- read.csv("Genes.csv")[2]
kept_genes <- read.csv("Kept Genes.csv")
kept_means <- read.csv("Kept Means.csv")
kept_standard_deviations <- read.csv("Kept Standard Deviations.csv")
gene_scores <- genedata_means + genedata_stdevs

