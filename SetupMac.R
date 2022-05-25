
library(tidyverse)
library(ggplot2)
genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_csv2 <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv")
drugdata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/primary-screen-replicate-collapsed-logfold-change.csv")
drugdata_csv2 <- t(data.frame(drugdata_csv)[-c(4,17,189,221,257,270,423,500,545,559),])
drugdata_matrix <- data.frame(matrix(as.numeric(matrix(drugdata_csv2, 4687, 568)[-c(1), ]), 4686, 568))
drugdata_drugs <- read.csv("Drugs.csv")[2]
drugdata_means <- read.csv("DrugMeans.csv")[2]
genedata_means <- read.csv("Means.csv")[2]
genedata_stdevs <- read.csv("Standard Deviations.csv")[2]
genedata_genes <- read.csv("Genes.csv")[2]
kept_genes <- read.csv("Kept Genes.csv")[2]
kept_means <- read.csv("Kept Means.csv")[2]
kept_drugs <- read.csv("Kept Drugs.csv")[2]
drugdata_sds <- read.csv("DrugSds.csv")[2]
kept_standard_deviations <- read.csv("Kept Standard Deviations.csv")[2]
gene_scores <- genedata_means + genedata_stdevs
kept_genes <- which(genedata_genes[,1] %in% kept_genes[,1])
kept_drugs <- which(drugdata_drugs[,1] %in% kept_drugs[,1])
genedata_csv3 <- genedata_csv2[,kept_genes]
drugdata_csv3 <- drugdata_csv2[c(1,kept_drugs),]

genedata_csv3[1,1]
genedata_cell_lines <- genedata_csv3[1,]
drugdata_cell_lines <- 