geneabridged <- function(genelist) {
  newlist <- genelist
  genelist1 <- genelist[,1]
  for (i in 1:length(genelist1)){
    newlist[i] <- substring(genelist1[i], 1, unlist(gregexpr('\\.', genelist1[i]))[1]-1)
  }
  return(newlist)
}

genedata_genes2 <- geneabridged(genedata_genes)


length(genedata_genes[,1])
substr(genedata_genes[,1][1], 1, unlist(gregexpr('\\.', genedata_genes[,1][1]))[1]-1)
length(genedata_genes[,1])