#This function removed the symbols in the genes that make them unreadable by the annotation tools.

geneabridged <- function(genelist) {
  newlist <- genelist
  genelist1 <- genelist
  for (i in 1:length(genelist1)){
    newlist[i] <- substring(genelist1[i], 1, unlist(gregexpr('\\.', genelist1[i]))[1]-1)
  }
  return(newlist)
}
#Testing if it worked
genedata_genes2 <- geneabridged(genedata_genes)
