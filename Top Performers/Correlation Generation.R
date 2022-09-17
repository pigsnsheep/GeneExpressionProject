count = 0


#Data formatting to make sure the for loop compiles correctly
  
genedrugnames <- rownames(correlations2)

classifications = annotations$drug_or_gene

unique(classifications)

correlations2 <- drop_na(data.frame(correlations))

correlations2$gene_or_drug <- classifications

classifications <- drop_na(data.frame(classifications))

rownames(genedrugannot) <- genedrugannot$V1

genedrugannot[genedrugnames[973],]$V4

#Initialize the vector

top_performers <- data.frame(matrix(ncol = 8, nrow = 0))

colnames(top_performers) <- c("drug", "drug_name","gene", "MoA", "panther_protein", "panther_family", "gene_of_interest", "correlation")


for (i in 1:3016) {
  for (j in 1:3016) {
    if(correlations[i,j]>0.5){
      if (classifications[i,] != classifications[j,]){
        if (classifications[i,] == "gene"){
          top_performers[nrow(top_performers)+1,] <- c(genedrugnames[j],  genedrugannot[genedrugnames[j],]$V4,  genedrugnames[i],  annotations$MoA[j],  annotations$panther_protein[i],  annotations$panther_family[i],  annotations$gene_of_interest[j], correlations[i,j])
        }
        else {
          top_performers[nrow(top_performers)+1,] <- c(genedrugnames[i],  genedrugannot[genedrugnames[i],]$V4,  genedrugnames[j],  annotations$MoA[i],  annotations$panther_protein[j],  annotations$panther_family[j],  annotations$gene_of_interest[i], correlations[i,j])
        }
      }
    }
  }
}

#Once I added the top performers, I manually removed pairs that didn't show covariation in the heatmap. I also ran this function with both the training and testing data



write.csv(top_performers, "top_performers.csv")
