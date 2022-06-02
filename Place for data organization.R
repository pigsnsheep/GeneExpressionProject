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
