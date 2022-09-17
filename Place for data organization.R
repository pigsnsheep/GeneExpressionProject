drugannot[,2]
kept_drugs
kept_drugs[,1]
for (i in 1:1078){
  kept_drugs[i,] <- substring(kept_drugs[i,],1, unlist(gregexpr('\\.',kept_drugs[i,]))[4]-1)
}
drugannotdrugs <- gsub('-','.',drugannot$broad_id)
write.csv(drugannotdrugs, "drugannot.csv")
drugannotdrugs <- which(kept_drugs %in% drugannotdrugs)
drugannotdrugs

rownames(genedatafinal) <- genedatafinal[,1]
genedatafinal <- genedatafinal[,2:559]
finalgenemeans <- rowMeans(genedatafinal)
finalgenes <- rownames(genedatafinal)
genedatafinal$sds <- rowSds(as.matrix(genedatafinal))
genedatafinal$Gene <- geneabridged(finalgenes)
genedatafinal2 <- genedatafinal %>% arrange(Gene, -sds) %>% filter(duplicated(Gene) == FALSE)
rownames(genedatafinal2) <- genedatafinal2$Gene
genedatafinal3 <- genedatafinal2[-c(559,560)]
write.csv(genedatafinal3, "Genedatafinal.csv")
