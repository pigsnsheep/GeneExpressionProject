#Generated cell line annotations

MakeDetermaIO <- function(testdata = matrix(), gmeans = NULL, gsds = NULL, logTransform = FALSE, scaled=T) {
  
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


MakeTNBCtype <- function(testsub = data.frame(),logTransform = F, scale=T, gmeans = NULL, gsds = NULL) {
  load("centroids.101.RData")
  centroids <- centroids.101[-which(rownames(centroids.101) %in% rownames(testsub) == F),]
  if(sum(rownames(testsub) %in% rownames(centroids) == F)>0) {testsub <- testsub[-which(rownames(testsub) %in% rownames(centroids) == F),]}
  
  #If needing to log transform
  if(logTransform){
    testsub[round(testsub,3)==0]<-0.001
    testsub<-log2(testsub )
  }
  
  centroids <- centroids[order(row.names(centroids)),]
  testsub <- testsub[order(row.names(testsub)),]
  
  # if given a gmeans, a gsds, or both, make sure they are formatted
  if(!is.null(gmeans)){ 
    gmeans <- merge(centroids,gmeans,by=0,all.x=TRUE)
    rownames(gmeans)<-gmeans[,1]
    gmeans<-gmeans[8]
    gmeans <- gmeans[order(rownames(gmeans)),,drop=FALSE]
    gmeans<-gmeans[,1]
  }
  if(!is.null(gsds)){
    gsds <- merge(centroids,gsds,by=0,all.x=TRUE)
    rownames(gsds)<-gsds[,1]
    gsds<-gsds[8]
    gsds <- gsds[order(rownames(gsds)),,drop=FALSE]
    gsds<-gsds[,1]
  }
  #perform only if you want to scale
  if (scale) {
    if(is.null(gmeans)){gmeans<-apply(testsub, 1, mean)}
    if(is.null(gsds)){gsds<-apply(testsub,1, sd)}
    testsub<-t(scale(t(testsub), gmeans, gsds))
  }
  #make sure all True
  
  headings101 <- c("Primary", "Secondary", "Hold",
                   "BL1_Corr", "BL2_Corr", "M_Corr","LAR_Corr", "MSL_Corr", "IM_101_Corr", 
                   "BL1_Cut", "BL2_Cut", "M_Cut","LAR_Cut", "MSL_Cut", "IM_101_Cut")
  TableRefSet <- t(testsub)
  TableRefSet <- TableRefSet[,1:15]
  TableRefSet[,1:2] <- "ND"
  TableRefSet [,3:15] <- 0.0000
  colnames(TableRefSet) <- headings101
  n <- ncol(testsub)
  UsingTable <- rep(0.000, n)  #variable to hold spearman data
  
  MatchStatement <- "No Match(es) at "
  for (i in 1:length(centroids)) {
    if(rownames(centroids)[i] != rownames(testsub)[i]) paste(MatchStatement, i, sep = " ")
  }
  if (MatchStatement == "No Match(es) at ") MatchStatement <- "Genes Match"
  
  print(MatchStatement)
  
  #perform Spearman 
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testsub[,i]), centroids$BL1, method = "spearman", "complete.obs")}
  TableRefSet[,4] <- UsingTable[1:n]
  
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testsub[,i]), centroids$BL2, method = "spearman", "complete.obs")}
  TableRefSet[,5] <- UsingTable[1:n]
  
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testsub[,i]), centroids$M, method = "spearman", "complete.obs")}
  TableRefSet[,6] <- UsingTable[1:n]
  
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testsub[,i]), centroids$LAR, method = "spearman", "complete.obs")}
  TableRefSet[,7] <- UsingTable[1:n]
  
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testsub[,i]), centroids$MSL, method = "spearman", "complete.obs")}
  TableRefSet[,8] <- UsingTable[1:n]
  
  for(i in 1:n) {UsingTable[i]<-cor(as.numeric(testsub[,i]), centroids$IM, method = "spearman", "complete.obs")}
  TableRefSet[,9] <- UsingTable[1:n]
  
  
  for(i in 4:9) {
    for (j in 1:nrow(TableRefSet)){
      if(as.numeric(TableRefSet[j,i]) > 0.19499){TableRefSet[j,i+6] <- 1}
    }
  }
  
  TableRefSet[,3] <- TableRefSet[,15]
  TableRefSet <- TableRefSet[,1:14]
  colnames(TableRefSet)[3] <- headings101[15]
  
  zDenom <-0.1428571 # sqrt of (1 / sqrt(101 -3))^2 + 1 / sqrt(101 -3))^2)
  
  for(i in 1:nrow(TableRefSet)) {
    SumAboveThreshold <- sum(as.numeric(TableRefSet[i,10:14]))
    if(SumAboveThreshold == 1) {
      locNames <- colnames(TableRefSet)[which.max(as.numeric(TableRefSet[i,4:8])) + 3]
      TableRefSet[i, 1] <- substr(locNames,1,regexpr("_", locNames)- 1)
    }
    
    if(SumAboveThreshold > 1) {
      #locNames <- colnames(TableRefSet)[which(rank(as.numeric(TableRefSet[i,4:8]))==5) + 3]
      #so I dither the data a bit
      set.seed(1)
      locNames <- colnames(TableRefSet)[which(rank(as.numeric(TableRefSet[i,4:8])+(runif(5)/10000))==5) + 3] # I changed this because sometimes there are ties
      locNames<-substr(locNames,1,regexpr("_", locNames)- 1)
      TableRefSet[i, 1] <-locNames 
      set.seed(1)
      pos1 <- as.numeric(which(rank(as.numeric(TableRefSet[i,4:8])+(runif(5)/10000)) ==5)) + 3
      set.seed(1)
      pos2 <- as.numeric(which(rank(as.numeric(TableRefSet[i,4:8])+(runif(5)/10000)) == 4)) +3
      
      GetZScore <- (as.numeric(TableRefSet[i,pos1]) -  as.numeric(TableRefSet[i,pos2])) / zDenom
      
      
      if (GetZScore < 1.96) {
        set.seed(1)
        secLocNames <- colnames(TableRefSet)[which(rank(as.numeric(TableRefSet[i,4:8])+(runif(5)/10000))==4) + 3]
        TableRefSet[i, 2] <- substr(secLocNames,1,regexpr("_", secLocNames)- 1)
      }
    }
  }
  
  return(TableRefSet)
  
}
