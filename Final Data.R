genedatafinal <- read.csv("Genedatafinal.csv", row.names=1)
cellannotationTNBC <- read.csv("cellannot.csv", row.names=1)
cellannotationDIO <- read.csv("DetermaIO-annot.csv", row.names=1)
geneannotation <- read.table("gene annotation.txt", sep="\t")


genedatafinal <- geneabridged(rownames(genedatafinal))
test <- MakeDetermaIO(genedatafinal3)
test2 <- MakeTNBCtype(genedatafinal3)
rownames(centroids.27)
order(rownames(centroids.27))
centroids <- centroids.27[order(rownames(centroids.27)),]
centroids <- centroids[-(which(row.names(centroids) %in%c("UBD", "KRT17", "IGJ"))),]

t(genedatafinal)[,1:2]

genedata_csv3 <- read.csv("Genedatafinal.csv")

test2 <- MakeTNBCtype(genedatafinal)

write.csv(test, "DetermaIO-annot.csv")
write.csv(test2, "cellannot.csv")

centroids <- centroids.101[which(rownames(centroids.101) %in% rownames(genedatafinal)),]

(rownames(centroids.101) %in% rownames(genedatafinal))
