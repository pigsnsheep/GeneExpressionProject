MakeDetermaIO <- function(testdata = matrix(), gmeans = NULL, gsds = NULL, logTransform = FALSE, scaled=TRUE) {
  
  centroids.27 <- load("centroids.27.RData")
  #centroids <- centroids.27[order(rownames(centroids.27)),]
  
  #if(length(rownames(centroids)) == 30) {centroids <- centroids[-(which(row.names(centroids) %in%c("UBD", "KRT17", "IGJ"))),]}
  
  testdata <- testdata[which(row.names(testdata) %in% row.names(centroids)),]
  centroids <- centroids[which(row.names(centroids) %in% row.names(testdata)),]
  
  centroids <- centroids[order(rownames(centroids)),]
  testdata <- testdata[order(rownames(testdata)),]
  
  if(!is.null(gmeans)){ # if given a gmeans, a gsds, or both, make sure they are formatted
    gmeans<-gmeans[rownames(gmeans) %in% rownames(centroids),,drop=FALSE]
    gmeans <- merge(centroids,gmeans,by=0,all.x=TRUE)
    rownames(gmeans)<-gmeans[,1]
    gmeans<-gmeans[4]
    gmeans <- gmeans[order(rownames(gmeans)),,drop=FALSE]
    
    gmeans<-gmeans[,1]
  }
  if(!is.null(gsds)){
    gsds<-gsds[rownames(gsds) %in% rownames(centroids),,drop=FALSE]
    gsds <- merge(centroids,gsds,by=0,all.x=TRUE)
    rownames(gsds)<-gsds[,1]
    gsds<-gsds[4]
    gsds <- gsds[order(rownames(gsds)),,drop=FALSE]
    gsds<-gsds[,1]
  }
  
  #If needing to log transform
  if(logTransform){
    testdata[round(testdata,3)==0]<-0.001
    testdata<-log2(testdata)
  }
  
  if (scaled) {
    if(is.null(gmeans)){gmeans<-apply(testdata, 1, mean)}
    if(is.null(gsds)){gsds<-apply(testdata,1, sd)}
    testdata<-t(scale(t(testdata), gmeans, gsds))
  }
  
  TableRefSet  <- NULL
  TableRefSet <- t(testdata)
  TableRefSet <- TableRefSet[,1:2]
  TableRefSet [,1:2] <- 0.0000
  colnames(TableRefSet) <- c("IM_Corr","IM_Cut")
  
  n <- ncol(testdata)
  
  UsingTable <- rep(0.000, n)
  
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testdata[,i]), centroids$IM, method = "spearman", "complete.obs")}
  TableRefSet[,1] <- UsingTable[1:n]
  
  TableRefSet[,2][which(TableRefSet[,1] >= 0.09)] <- 1
  
  return(TableRefSet)
}
