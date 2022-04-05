
library(tidyverse)
library(ggplot2)
genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_csv2 <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv")
drugdata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/primary-screen-replicate-collapsed-logfold-change.csv")
drugdata_csv2 <- t(data.frame(drugdata_csv)[-c(4,17,189,221,257,270,423,500,545,559),])
drugdata_matrix <- data.frame(matrix(as.numeric(matrix(drugdata_csv2, 4687, 568)[-c(1), ]), 4686, 568))
drugdata_drugs <- rownames(drugdata_csv2)
drugdata_means <- rowMeans(drugdata_matrix, na.rm = TRUE)
genedata_means <- read.csv("Means.csv")[2]
genedata_stdevs <- read.csv("Standard Deviations.csv")[2]
genedata_genes <- read.csv("Genes.csv")[2]
kept_genes <- read.csv("Kept Genes.csv")[2]
kept_means <- read.csv("Kept Means.csv")[2]
kept_standard_deviations <- read.csv("Kept Standard Deviations.csv")[2]
gene_scores <- genedata_means + genedata_stdevs

sum(duplicated(drugdata_drugs) == TRUE)
mean(drugdata_csv2[2,], na.rm = TRUE)
drugdata_csv2[,1]

write.csv(drugdata_drugs, "Drugs.csv")
drugdata_sds <- apply(drugdata_matrix, 1, sd, na.rm = TRUE)
write.csv(drugdata_means, "DrugMeans.csv")
write.csv(drugdata_sds, "DrugSds.csv")
