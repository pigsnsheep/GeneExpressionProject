
library(tidyverse)
library(ggplot2)
genedata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv", header = TRUE)
genedata_csv2 <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/CCLE_expression.csv")
drugdata_csv <- read.csv("/Users/pigsnsheep/Downloads/Gene\ Expression\ Project/primary-screen-replicate-collapsed-logfold-change.csv")
drugannot <- na.omit(read.csv("primary-screen-replicate-treatment-info.csv"))
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
genedata_genes1 <- geneabridged(genedata_genes)
kept_gene_indices <- which(unique(geneabridged(genedata_genes)) %in% kept_genes[,1])
kept_drug_indices <- which(drugdata_drugs[,1] %in% kept_drugs[,1])
genedata_csv3 <- genedata_csv2[,2:19178][,kept_gene_indices]
drugdata_csv3 <- drugdata_csv2[c(1,kept_drugs),]
common_cell_lines <- genedata_cell_lines[which(drugdata_cell_lines %in% genedata_cell_lines)]
common_cell_lines_1 <- which(drugdata_cell_lines %in% genedata_cell_lines)
common_cell_lines_2 <- which(genedata_cell_lines %in% drugdata_cell_lines)
genedata_csv4 <- genedata_csv3[common_cell_lines_2,]
drugdata_csv4 <- drugdata_csv3[,common_cell_lines_1]
colnames(drugdata_csv3) <- colnames(genedata_csv3)
compound_data <- rbind(genedata_csv4, drugdata_csv4)

genedata_cell_lines <- genedata_csv2[,1]
drugdata_cell_lines <- drugdata_csv2[1,]

gsub("\\.",",",genedata_genes[1,1])

rownames(genedata_csv4) <- kept_genes[,1]
length(kept_genes)

kept_genes_3 <- kept_genes[kept_gene_indices,]

rownames(genedata_csv4) <- common_cell_lines

genedata_csv5 <- t(genedata_csv4)



write.csv(genedata_csv5, "Genedatafinal.csv")


