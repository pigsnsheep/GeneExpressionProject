genedatafinal <- read.csv("Genedatafinal.csv")
rownames(genedatafinal) <- genedatafinal[,1]
genedatafinal <- genedatafinal[,2:559]
finalgenemeans <- rowMeans(genedatafinal)

test <- MakeDetermaIO(genedatafinal)
rownames(centroids.27)
order(rownames(centroids.27))
centroids <- centroids.27[order(rownames(centroids.27)),]
centroids <- centroids[-(which(row.names(centroids) %in%c("UBD", "KRT17", "IGJ"))),]

t(genedatafinal)[,1:2]
