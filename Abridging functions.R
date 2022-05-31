geneabridged <- function(genelist) {
  for (i in 1:length(genelist)){
    genelist[i] <- substr(genelist[i], 1, unlist(gregexpr('\\.', genelist[i]))[1])
  }
}