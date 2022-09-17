cluster_plus_ultra<-function(mapData=NULL, 
                             tissueannot=NULL, #the first row is 0 or 1, depending on whether to use that columns annotation
                             geneannot=NULL, #the first row is 0 or 1, depending on whether to use that columns annotation
                             scale.genes=TRUE, scale.samples=TRUE,
                             annotate_genes=FALSE,
                             cluster.method="sumsquares", col_cluster_method="kmeans", ordering_columns="diagonal", 
                             return.plot=TRUE, return.stats=FALSE, return_data=FALSE) {
  ## mapData is a matrix of genes (rows) and tissues (columns)
  ## tissueannot and geneannot contain annotation information. The first row is 0 or 1, and determines whether the column is used
  ## cluster.method is either kmeans, 'pearson', or 'sumsquares' (default)
  ## col_cluster_method is either kmeans (default) or 101centroid
  ## ordering_columns is either 'right' (the original method) or "diagonal" (the default)
  ## return.plot determines whether to return the plot
  ## return.stats determines whether to return some summary stats
  ## return_data determines whether to return the clustered data
  
  require(ComplexHeatmap)
  require(RColorBrewer)
  require(circlize) #for any annotations uses colorRamp2
  
  holdRowNames <- rownames(mapData)
  
  mapData<-mapData[,order(colnames(mapData))]  #lets order things to be safe
  #create the subtypevector from the tissue annotation
  SubtypeVector<-tissueannot[-1,]
  colnames(SubtypeVector)<-tissueannot[1,]
  SubtypeVector<-SubtypeVector[-1,]
  SubtypeVector<-SubtypeVector[,colnames(SubtypeVector) %in% c("IM_Corr","MSL_Corr","M_Corr")]
  SubtypeVector<-SubtypeVector[rownames(SubtypeVector) %in% colnames(mapData),]
  SubtypeVector<-SubtypeVector[order(rownames(SubtypeVector)),]
  therownames<-rownames(SubtypeVector)
  SubtypeVector$IM_Corr <- lapply(SubtypeVector$IM_Corr, remove_first_decimal)
  SubtypeVector <- apply(SubtypeVector, 2,function(x) as.numeric(as.character(x)))
  rownames(SubtypeVector)<-therownames
  
  unscaledSTV<-SubtypeVector # save the unscaled vector in case its needed.
  SubtypeVector <- scale(SubtypeVector)
  if (scale.samples) {mapData <- scale(mapData)}
  #mapData <- rbind(mapData, t(SubtypeVector))
  if (scale.genes) { mapData <- t(scale(t(mapData)))}
  
  SubtypeVector <- t(mapData[(nrow(mapData)-2):nrow(mapData),])
  #
  #mapData <- mapData[1:(nrow(mapData)-3), ]
  
  mapData[is.na(mapData)] <- 0
  
  mdMappingData <- list(mappingData = matrix(), 
                        row.size = c(0,0,0), 
                        col.size = c(0,0,0), 
                        col.position  = c(left = 1, middle = 2, right = 3))
  
  rCl<- integer(nrow(mapData)) #create a named vector
  names(rCl)<-rownames(mapData)
  cCl<- integer(ncol(mapData)) #create a named vector
  names(cCl)<-colnames(mapData)
  kRows<- NULL
  kColumns<-NULL
  
  if(cluster.method=="pearson") {
    for (i in 1:nrow(mapData)) {
      ds<-apply(SubtypeVector, 2, function(x) cor(x, t(mapData)[,i],method = "pearson"))
      rCl[i]<-match(max(ds),ds)
    }
    kRows$size[1]<-sum(rCl == 1)
    kRows$size[2]<-sum(rCl == 2)
    kRows$size[3]<-sum(rCl == 3)
  }
  
  if(cluster.method=="sumsquares") {
    
    for (i in 1:nrow(mapData)) {
      ds<-apply(SubtypeVector, 2, function(x) sum((t(mapData)[,i]-x)^2))
      rCl[i]<-match(min(ds),ds)
    }
    
    kRows$size[1]<-sum(rCl == 1)
    kRows$size[2]<-sum(rCl == 2)
    kRows$size[3]<-sum(rCl == 3)
  }
  
  if(cluster.method=="kmeans") {
    set.seed(2)
    kRows <- kmeans(mapData, centers = t(SubtypeVector[,1:3]))
    rCl <- kRows$cluster
  }
  
  if(col_cluster_method=="101centroid") {
    for (i in 1:ncol(mapData)) {
      ds<-SubtypeVector[colnames(mapData)[i],]
      #ds<-unscaledSTV[colnames(mapData)[i],] # if you want the unscaled vector
      cCl[i]<-match(max(ds),ds)
    }
    kColumns$size[1]<-sum(cCl == 1)
    kColumns$size[2]<-sum(cCl == 2)
    kColumns$size[3]<-sum(cCl == 3)
  }
  
  if(col_cluster_method=="kmeans") {
    set.seed(2)
    #set.seed(NULL)
    
    kColumns <- kmeans(t(mapData), centers = 3)
    cCl <- kColumns$cluster
  }
  
  #a cludgy fix for clusters with zero members
  if(kRows$size[1]==0) {
    if(kRows$size[2]>3){kRows$size<-c(2,kRows$size[2]-2,kRows$size[3])} else {
      kRows$size<-c(2,kRows$size[2],kRows$size[3]-2)
    }
  }
  if(kRows$size[2]==0) {
    if(kRows$size[1]>2){kRows$size<-c(kRows$size[1]-2,2,kRows$size[3])} else {
      kRows$size<-c(kRows$size[1],2,kRows$size[3]-2)
    }
  }
  if(kRows$size[3]==0) {
    if(kRows$size[1]>2){kRows$size<-c(kRows$size[1]-2,kRows$size[2],2)} else {
      kRows$size<-c(kRows$size[1],kRows$size[2]-2,2)
    }
  }
  
  mdMappingData$row.size <- kRows$size
  mdMappingData$col.size <- kColumns$size
  
  md2 <- mapData[c(which(rCl == 1), which(rCl == 2), which(rCl == 3)),  
                 c(which(cCl == 1), which(cCl == 2), which(cCl == 3))]
  
  nRowSect_1 <- kRows$size[1]  #number of rows in IM
  nRowSect_2 <- kRows$size[2] #number of rows in MSL
  nRowSect_3 <- kRows$size[3]  #number of rows in M
  
  rSectCoor_1 <- c(1,kRows$size[1])
  rSectCoor_2 <- c(rSectCoor_1[2] + 1,rSectCoor_1[2] + kRows$size[2])
  rSectCoor_3 <- c(rSectCoor_2[2] + 1,rSectCoor_2[2] + kRows$size[3])
  
  mdRowSect_1<- md2[c(hclust(dist(md2[rSectCoor_1[1]: rSectCoor_1[2],]))$order),]  #cluster IM within IM
  mdRowSect_2<- md2[rSectCoor_1[2] + c(hclust(dist(md2[rSectCoor_2[1]:rSectCoor_2[2],]))$order),]  #MSL
  mdRowSect_3<- md2[rSectCoor_2[2] + c(hclust(dist(md2[rSectCoor_3[1]:rSectCoor_3[2],]))$order),] #M
  
  #if no clustering within clusters-THIS has issues
  # mdRowSect_1<- md2[rSectCoor_1[1]: rSectCoor_1[2],]  #cluster IM
  # mdRowSect_2<- md2[rSectCoor_2[1]:rSectCoor_2[2],] #MSL
  # mdRowSect_3<- md2[rSectCoor_3[1]:rSectCoor_3[2],] #M
  
  md2 <- rbind(mdRowSect_1, mdRowSect_2, mdRowSect_3)
  md2 <- t(md2)
  
  nColSect_1 <- kColumns$size[1]  #number of columns in k = 1
  nColSect_2 <- kColumns$size[2]  #number of columns in k = 2
  nColSect_3 <- kColumns$size[3]  #number of columns in k = 3
  
  cSectCoor_1 <- c(1, kColumns$size[1])
  cSectCoor_2 <- c(cSectCoor_1[2] + 1, cSectCoor_1[2] + kColumns$size[2])
  cSectCoor_3 <- c(cSectCoor_2[2] + 1, cSectCoor_2[2] + kColumns$size[3])
  
  mdColSect_1<- md2[c(hclust(dist(md2[cSectCoor_1[1]:cSectCoor_1[2],]))$order),]  #cluster 1 within 1
  mdColSect_2<- md2[nColSect_1 + 
                      c(hclust(dist(md2[cSectCoor_2[1]:cSectCoor_2[2],]))$order),]  #2
  mdColSect_3<- md2[nColSect_1 +nColSect_2 + 
                      c(hclust(dist(md2[cSectCoor_3[1]:cSectCoor_3[2],]))$order),] #3
  
  ## if no clustering within clusters:
  # mdColSect_1<- md2[cSectCoor_1[1]:cSectCoor_1[2],]  #cluster 1
  # mdColSect_2<- md2[cSectCoor_2[1]:cSectCoor_2[2],]  #2
  # mdColSect_3<- md2[cSectCoor_3[1]:cSectCoor_3[2],] #3
  # # 
  md2 <- rbind(mdColSect_1, mdColSect_2, mdColSect_3)
  mdColSectComb <- list(mdColSect_1, mdColSect_2 ,mdColSect_3)
  
  ## order the three columns so they highlight the TIME structure
  if(ordering_columns=="diagonal"){
    require(gtools)
    arraylist<-permutations(3,3,c(1:3))
    expsum<-NULL
    for (i in 1:nrow(arraylist)) {
      expsum[i]<-eval(parse(text=(paste0("mean(mdColSect_",arraylist[i,1],"[,rSectCoor_1[1]:rSectCoor_1[2]]) + mean(mdColSect_",arraylist[i,2],"[,rSectCoor_2[1]:rSectCoor_2[2]]) + mean(mdColSect_",arraylist[i,3],"[,rSectCoor_3[1]:rSectCoor_3[2]])"))))
    }
    left<-arraylist[order(expsum,decreasing = TRUE)[1],1]
    middle<-arraylist[order(expsum,decreasing = TRUE)[1],2]
    right<-arraylist[order(expsum,decreasing = TRUE)[1],3]
  }
  
  if(ordering_columns=="right"){
    left <- which(rank(c(mean(mdColSect_1[,rSectCoor_1[1]:rSectCoor_1[2]]), 
                         mean(mdColSect_2[,rSectCoor_1[1]:rSectCoor_1[2]]) ,
                         mean(mdColSect_3[,rSectCoor_1[1]:rSectCoor_1[2]]))) 
                  == 3) 
    right <- which(rank(c(mean(mdColSect_1[,rSectCoor_3[1]:rSectCoor_3[2]]),
                          mean(mdColSect_2[,rSectCoor_3[1]:rSectCoor_3[2]]) ,
                          mean(mdColSect_3[,rSectCoor_3[1]:rSectCoor_3[2]])))
                   == 3)
    if(left == right) {
      #should not happen so print a warning
      print('Warning: Highest mean for IM and M are same cluster')
      right <- which(rank(c(mean(mdColSect_1[,rSectCoor_3[1]:rSectCoor_3[2]]),
                            mean(mdColSect_2[,rSectCoor_3[1]:rSectCoor_3[2]]) ,
                            mean(mdColSect_3[,rSectCoor_3[1]:rSectCoor_3[2]])))
                     == 1)
    }
    middle <- setdiff(c(1:3), c(left,right))  #remaining group goes in middle
  }
  
  #reorder columns based on this ordering
  md2 <- rbind(data.frame(mdColSectComb[left]),data.frame(mdColSectComb[middle]),data.frame(mdColSectComb[right]))
  
  nColSect_1 <- length(which(cCl == left))  #number of columns in k = left most sector
  nColSect_2 <- length(which(cCl == middle))  #number of columns in k = middle sector
  nColSect_3 <- length(which(cCl == right))  #number of columns in k = right most sector
  
  cSectCoor_1 <- c(1,nColSect_1)
  cSectCoor_2 <- c(nColSect_1 + 1,nColSect_1 + nColSect_2)
  cSectCoor_3 <- c(nColSect_1 + nColSect_2 + 1,nrow(md2))
  
  mdMappingData$col.position["left"] <- left
  mdMappingData$col.position["middle"] <- middle
  mdMappingData$col.position["right"] <- right
  
  #prep md2 for heatmap
  md2 <- t(md2)
  mdMappingData$mappingData <- md2
  class(mdMappingData) <- append(class(mdMappingData), "make_md2")
  md2<-mdMappingData
  
  #### annotations
  ################
  
  ### row (gene) annotations
  #############
  ## create vector of gene symbols in data to be clustered
  rowsincluster<-data.frame(row.names(md2$mappingData))
  colnames(rowsincluster)<-"clusterID"
  rowsincluster$order_id  <- 1:nrow(rowsincluster)
  
  rHa1<-NULL
  if(annotate_genes) {
    ## create the annotation list
    #the first row of geneannot determines whether or not to use the column
    geneannot<-geneannot[,!(geneannot[1,]==0)]
    geneannot<-geneannot[-1,]
    colnames(geneannot)<-geneannot[1,]
    geneannot<-geneannot[-1,]
    geneannot <-merge(x = geneannot, y = rowsincluster, by.x = "gene", by.y="clusterID", all.y=TRUE)
    geneannot<-geneannot[order(geneannot$order_id),]
    
    #convert appropriate columns to numeric
    for(i in 1:ncol(geneannot)){
      if(sum(!grepl('^-?[0-9.]+$', na.omit(geneannot[,i])))==0){
        geneannot[,i]<-as.numeric(as.numeric(geneannot[,i]))
      }
    }
    
    annotcolors<-vector(mode = "list", length = ncol(geneannot)-1)
    rhastatement<-"rHa1<- rowAnnotation("
    colstatement<-"col=list("
    
    for(annotation_count in 1:(ncol(geneannot)-2)){ #don't grab the last column, which is used for ordering
      #don't grab the last column, which is used for ordering
      #identify top set to limit annotation, if needed
      if(!is.numeric(geneannot[,annotation_count+1]) && length(unique(geneannot[,annotation_count+1]))>25){
        topannots<-data.frame(sort(table(geneannot[,annotation_count+1]),decreasing=TRUE)[1:25])
        is.na(geneannot[,annotation_count+1]) <- !(geneannot[,annotation_count+1] %in% topannots$Var1)
      }
      #then make the colors
      
      #if numeric and a lot of values, then use a color ramp of red>yellow>blue
      if(is.numeric(geneannot[,annotation_count+1]) && 
         length(unique(na.omit(geneannot[,annotation_count+1])))>20){
        tempannot<-geneannot[,annotation_count+1]
        
        #if the range is huge, limit it prior to determining colors
        colorsdfilter<-3
        if(min(tempannot,na.rm = TRUE)<(mean(tempannot,na.rm = TRUE)-colorsdfilter*sd(tempannot,na.rm = TRUE)) || max(tempannot,na.rm = TRUE)>(mean(tempannot,na.rm = TRUE)-colorsdfilter*sd(tempannot,na.rm = TRUE))) {
          #tempannot<-log2(tempannot-min(tempannot,na.rm=TRUE)+1)
          tempannot<-tempannot[tempannot<(mean(tempannot,na.rm = TRUE)+colorsdfilter*sd(tempannot,na.rm = TRUE)) & tempannot>(mean(tempannot,na.rm = TRUE)-colorsdfilter*sd(tempannot,na.rm = TRUE))]
        }
        annotcolors[[annotation_count]] = colorRamp2(c(min(tempannot,na.rm = TRUE), (min(tempannot,na.rm = TRUE)+max(tempannot,na.rm = TRUE))/2, max(tempannot,na.rm = TRUE)), c("red", "yellow", "blue"))
      }
      # if all NA's, then make white
      if(sum(is.na(geneannot[,annotation_count+1]))==length(geneannot[,annotation_count+1])){
        geneannot[,annotation_count+1]<-999
        annotcolors[[annotation_count]] <-  c('999' = "white")
      }
      #if numeric and not too many values, use a color ramp of white to dark blue
      if(is.numeric(geneannot[,annotation_count+1]) && length(unique(na.omit(geneannot[,annotation_count+1])))<=20 && length(unique(na.omit(geneannot[,annotation_count+1])))>1){
        annotcolors[[annotation_count]] = colorRamp2(c(min(geneannot[,annotation_count+1],na.rm = TRUE), max(geneannot[,annotation_count+1],na.rm = TRUE)), c("white","darkblue"))
      }
      
      #if has the values 0 & 1, or L & H, then choose two colors
      if(sum(!grepl('^-?[01LH]+$', na.omit(geneannot[,annotation_count+1])))==0){
        if(length(unique(na.omit(geneannot[,annotation_count+1])))<=2){
          if('L' %in% geneannot[,annotation_count+1] || 'H' %in% geneannot[,annotation_count+1]) {
            annotcolors[[annotation_count]] <-  c('L' = "lightgoldenrod", 'H' = "skyblue")} 
          if('0' %in% geneannot[,annotation_count+1] || '1' %in% geneannot[,annotation_count+1]) {
            annotcolors[[annotation_count]] <-  c('0' = "lightgoldenrod", '1' = "skyblue")} 
        }
      }
      
      #if non-numeric, then use colorRampPalette
      if(!is.numeric(geneannot[,annotation_count+1])){
        annotcolors[[annotation_count]] <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(factor(geneannot[,annotation_count+1]))))
        names(annotcolors[[annotation_count]]) <- levels(factor(geneannot[,annotation_count+1]))
      }
      
      #decide whether its a factor or not
      if(!is.numeric(geneannot[,annotation_count+1])) {
        rhastatement<-paste(rhastatement,paste0(colnames(geneannot)[annotation_count+1]," = factor(geneannot[,",annotation_count+1,"]),"),sep=" ")}
      if(is.numeric(geneannot[,annotation_count+1])) {
        rhastatement<-paste(rhastatement,paste0(colnames(geneannot)[annotation_count+1]," = geneannot[,",annotation_count+1,"],"),sep=" ")}
      
      colstatement<-paste(colstatement,paste0(colnames(geneannot)[annotation_count+1]," = annotcolors[[",annotation_count,"]],"),sep=" ")
    }
    colstatement<-substr(colstatement,1,nchar(colstatement)-1) #remove the last comma from colstatement
    colstatement<-paste0(colstatement,"), na_col = 'white')")
    rhastatement<-paste(rhastatement,colstatement,sep=" ")
    eval(parse(text=rhastatement))
  }
  ############## 
  
  ### column (sample) annotations
  #############
  
  mdtemp<-as.data.frame(rownames(t(md2[["mappingData"]])))
  colnames(mdtemp)<-"sampleid"
  mdtemp$order_id<-1:nrow(mdtemp) #make sure we have the order of the columns after clustering
  
  #first row of annotation determines which columns to use
  tissueannot<-tissueannot[,!(tissueannot[1,]==0)]
  tissueannot<-tissueannot[-1,]
  colnames(tissueannot)<-tissueannot[1,]
  tissueannot<-tissueannot[-1,]
  write.csv(tissueannot, "test123.csv")
  write.csv(mdtemp, 'test1234.csv')
  tissueannot$sampleid <- rownames(tissueannot)
  tissueannot <-merge(x = tissueannot, y = mdtemp, by = "sampleid", all.y=TRUE)
  tissueannot<-tissueannot[order(tissueannot$order_id),]
  tissueannot<-as.data.frame(tissueannot)
  
  #convert appropriate columns to numeric
  for(i in 1:ncol(tissueannot)){
    if(sum(!grepl('^-?[0-9.]+$', na.omit(tissueannot[,i])))==0){
      tissueannot[,i]<-as.numeric(as.numeric(tissueannot[,i]))
    }
  }
  
  annotcolors<-vector(mode = "list", length = ncol(tissueannot)-1)
  chastatement<-"colAnn <- HeatmapAnnotation( which = 'col',"
  colstatement<-"col=list(type = c("
  
#  for(annotation_count in 1:(ncol(tissueannot)-2)){ #don't grab the last column, which is used for ordering
  for(annotation_count in 1:2){ #don't grab the last column, which is used for ordering
  #identify top set to limit annotation, if needed
    if(!is.numeric(tissueannot[,annotation_count+1]) && length(unique(tissueannot[,annotation_count+1]))>35){
      topannots<-data.frame(sort(table(tissueannot[,annotation_count+1]),decreasing=TRUE)[1:35])
      is.na(tissueannot[,annotation_count+1]) <- !(tissueannot[,annotation_count+1] %in% topannots$Var1)
    }
    #then make the colors
    
    #if numeric and a lot of values, then use a color ramp of red>yellow>blue
    if(is.numeric(tissueannot[,annotation_count+1]) && 
       length(unique(na.omit(tissueannot[,annotation_count+1])))>20){
      tempannot<-tissueannot[,annotation_count+1]
      
      #if the range is huge, limit it prior to determining colors
      colorsdfilter<-3
      if(min(tempannot,na.rm = TRUE)<(mean(tempannot,na.rm = TRUE)-colorsdfilter*sd(tempannot,na.rm = TRUE)) || max(tempannot,na.rm = TRUE)>(mean(tempannot,na.rm = TRUE)-colorsdfilter*sd(tempannot,na.rm = TRUE))) {
        #tempannot<-log2(tempannot-min(tempannot,na.rm=TRUE)+1)
        tempannot<-tempannot[tempannot<(mean(tempannot,na.rm = TRUE)+colorsdfilter*sd(tempannot,na.rm = TRUE)) & tempannot>(mean(tempannot,na.rm = TRUE)-colorsdfilter*sd(tempannot,na.rm = TRUE))]
      }
      annotcolors[[annotation_count]] = colorRamp2(c(min(tempannot,na.rm = TRUE), (min(tempannot,na.rm = TRUE)+max(tempannot,na.rm = TRUE))/2, max(tempannot,na.rm = TRUE)), c("red", "yellow", "blue"))
    }
    # if all NA's, then make white
    if(sum(is.na(tissueannot[,annotation_count+1]))==length(tissueannot[,annotation_count+1])){
      tissueannot[,annotation_count+1]<-999
      annotcolors[[annotation_count]] <-  c('999' = "white")
    }
    #if numeric and not too many values, use a color ramp of white to dark blue
    if(is.numeric(tissueannot[,annotation_count+1]) && length(unique(na.omit(tissueannot[,annotation_count+1])))<=20 && length(unique(na.omit(tissueannot[,annotation_count+1])))>1){
      annotcolors[[annotation_count]] = colorRamp2(c(min(tissueannot[,annotation_count+1],na.rm = TRUE), max(tissueannot[,annotation_count+1],na.rm = TRUE)), c("white","darkblue"))
    }
    
    #if has the values 0 & 1, or L & H, then choose two colors
    if(sum(!grepl('^-?[01LH]+$', na.omit(tissueannot[,annotation_count+1])))==0){
      if(length(unique(na.omit(tissueannot[,annotation_count+1])))<=2){
        if('L' %in% tissueannot[,annotation_count+1] || 'H' %in% tissueannot[,annotation_count+1]) {
          annotcolors[[annotation_count]] <-  c('L' = "lightgoldenrod", 'H' = "skyblue")} 
        if('0' %in% tissueannot[,annotation_count+1] || '1' %in% tissueannot[,annotation_count+1]) {
          annotcolors[[annotation_count]] <-  c('0' = "lightgoldenrod", '1' = "skyblue")} 
      }
    }
    
    #if non-numeric, then use colorRampPalette
    if(!is.numeric(tissueannot[,annotation_count+1])){
      annotcolors[[annotation_count]] <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(factor(tissueannot[,annotation_count+1]))))
      names(annotcolors[[annotation_count]]) <- levels(factor(tissueannot[,annotation_count+1]))
    }
    
    #decide whether its a factor or not
    if(!is.numeric(tissueannot[,annotation_count+1])) {
      colstatement<-paste(colstatement,paste0(colnames(tissueannot)[annotation_count+1]," = factor(tissueannot[,",annotation_count+1,"]),"),sep=" ")}
    if(is.numeric(tissueannot[,annotation_count+1])) {
      colstatement<-paste(colstatement,paste0("tissueannot[,",annotation_count+1,"] = ", colnames(tissueannot)[annotation_count+1],", "),sep=" ")}
    
    colstatement<-paste(colstatement,paste0("annotcolors[[",annotation_count,"]] = ", colnames(tissueannot)[annotation_count+1], ", "),sep=" ")
  }
  colstatement<-substr(colstatement,1,nchar(colstatement)-1) #remove the last comma from colstatement
  colstatement<-paste0(colstatement,"), na_col = 'white')")
  chastatement<-paste(chastatement,colstatement,")",sep=" ")
  print(chastatement)
  
  print(md2)
  print(tissueannot)
  
  #XXX: Testing eval(parse(text=chastatement))
  colAnn <- HeatmapAnnotation( 
    df = md2$mappingData, 
    na_col = 'white')
  
  #and now lets create the plot object, if called for
  clusterplot<-NULL
  if (return.plot){
    TisMap <- suppressMessages(Heatmap(
      md2$mappingData,
      name = "a trial run", 
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 6),
      show_column_names = FALSE,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      top_annotation= NULL,
      #top_annotation=colAnn,
      right_annotation = rHa1,
      #right_annotation = NULL,
      row_split = factor(c(rep("IM", md2$row.size[1]),
                           rep("MSL", md2$row.size[2]),
                           rep("M", md2$row.size[3])),
                         levels= c("IM","MSL", "M" )),
      
      column_split = factor(c(rep("(1)", md2$col.size[md2$col.position["left"]]),
                              rep("(2)", md2$col.size[md2$col.position["middle"]]),
                              rep("(3)", md2$col.size[md2$col.position["right"]])), 
                            levels= c("(1)", "(2)", "(3)")),
      
      column_title_side = "bottom",
      show_column_dend = F,
      cluster_column_slices = FALSE,
      cluster_row_slices = FALSE,
      show_row_dend = FALSE,
      show_heatmap_legend = F
    ))
    
    clusterplot<-TisMap
  }
  
  ## lets add in the clusters to the gene and tissue annotation before returning it
  geneannot$cluster<-c(rep("IM",md2$row.size[1]),rep("MSL",md2$row.size[2]),rep("M",md2$row.size[3]))
  tissueannot$cluster<-c(rep("1",md2$col.size[md2$col.position["left"]]),rep("2",md2$col.size[md2$col.position["middle"]]),rep("3",md2$col.size[md2$col.position["right"]]))
  
  #return the cluster data, if called for
  clustered_data<-NULL
  if(return_data) {clustered_data<-md2$mappingData}
  
  #the next section calculates some statistics, if called for
  IM.prop<-NULL
  gene.dist<-NULL
  sample.dist<-NULL
  row.pos<-NULL
  col.pos<-NULL
  cluster.avg<-NULL
  tissue_part_summary<-NULL
  
  col.pos<-c(md2$col.size[md2$col.position["left"]], 
             md2$col.size[md2$col.position["middle"]],
             md2$col.size[md2$col.position["right"]])
  row.pos<-md2$row.size
  
  if(return.stats){
    ## and some summary data
    ## lets determine distance of the clusters from each other.
    ## this kind of got long. All the vectors have different lengths, 
    ## so it's tricky to write this well. I bet there was a cleverer way.
    ## The cluster is represented as 9 clusters labeled 1 through 9, 1=IM-1, 2=IM-2, 4=MSL-1, etc.
    
    gene.centroid.1<-colMeans(md2$mappingData[1:row.pos[1],1:col.pos[1]])
    gene.centroid.4<-colMeans(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),1:col.pos[1]])
    gene.centroid.7<-colMeans(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),1:col.pos[1]])
    
    gene.centroid.2<-colMeans(md2$mappingData[1:row.pos[1],col.pos[1]:(col.pos[1]+col.pos[2])])
    gene.centroid.5<-colMeans(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),col.pos[1]:(col.pos[1]+col.pos[2])])
    gene.centroid.8<-colMeans(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),col.pos[1]:(col.pos[1]+col.pos[2])])
    
    gene.centroid.3<-colMeans(md2$mappingData[1:row.pos[1],(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    gene.centroid.6<-colMeans(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    gene.centroid.9<-colMeans(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    
    sample.centroid.1<-rowMeans(md2$mappingData[1:row.pos[1],1:col.pos[1]])
    sample.centroid.4<-rowMeans(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),1:col.pos[1]])
    sample.centroid.7<-rowMeans(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),1:col.pos[1]])
    
    sample.centroid.2<-rowMeans(md2$mappingData[1:row.pos[1],col.pos[1]:(col.pos[1]+col.pos[2])])
    sample.centroid.5<-rowMeans(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),col.pos[1]:(col.pos[1]+col.pos[2])])
    sample.centroid.8<-rowMeans(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),col.pos[1]:(col.pos[1]+col.pos[2])])
    
    sample.centroid.3<-rowMeans(md2$mappingData[1:row.pos[1],(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    sample.centroid.6<-rowMeans(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    sample.centroid.9<-rowMeans(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    
    gene.centroid.list<-list(gene.centroid.1,gene.centroid.2,gene.centroid.3,gene.centroid.4,gene.centroid.5,gene.centroid.6,gene.centroid.7,gene.centroid.8,gene.centroid.9)
    sample.centroid.list<-list(sample.centroid.1,sample.centroid.2,sample.centroid.3,sample.centroid.4,sample.centroid.5,sample.centroid.6,sample.centroid.7,sample.centroid.8,sample.centroid.9)
    
    gene.dist<-array(NA,dim=c(9,3))
    gene.dist[1,]<-c(dist(rbind(gene.centroid.list[[1]],gene.centroid.list[[4]])),1,4)
    gene.dist[2,]<-c(dist(rbind(gene.centroid.list[[1]],gene.centroid.list[[7]])),1,7)
    gene.dist[3,]<-c(dist(rbind(gene.centroid.list[[4]],gene.centroid.list[[7]])),4,7)
    gene.dist[4,]<-c(dist(rbind(gene.centroid.list[[2]],gene.centroid.list[[5]])),2,5)
    gene.dist[5,]<-c(dist(rbind(gene.centroid.list[[2]],gene.centroid.list[[8]])),2,8)
    gene.dist[6,]<-c(dist(rbind(gene.centroid.list[[5]],gene.centroid.list[[8]])),5,8)
    gene.dist[7,]<-c(dist(rbind(gene.centroid.list[[3]],gene.centroid.list[[6]])),3,6)
    gene.dist[8,]<-c(dist(rbind(gene.centroid.list[[3]],gene.centroid.list[[9]])),3,9)
    gene.dist[9,]<-c(dist(rbind(gene.centroid.list[[6]],gene.centroid.list[[9]])),6,9)
    
    sample.dist<-array(NA,dim=c(9,3))
    sample.dist[1,]<-c(dist(rbind(sample.centroid.list[[1]],sample.centroid.list[[2]])),1,2)
    sample.dist[2,]<-c(dist(rbind(sample.centroid.list[[1]],sample.centroid.list[[3]])),1,3)
    sample.dist[3,]<-c(dist(rbind(sample.centroid.list[[2]],sample.centroid.list[[3]])),2,3)
    sample.dist[4,]<-c(dist(rbind(sample.centroid.list[[4]],sample.centroid.list[[5]])),4,5)
    sample.dist[5,]<-c(dist(rbind(sample.centroid.list[[4]],sample.centroid.list[[6]])),4,6)
    sample.dist[6,]<-c(dist(rbind(sample.centroid.list[[5]],sample.centroid.list[[6]])),5,6)
    sample.dist[7,]<-c(dist(rbind(sample.centroid.list[[7]],sample.centroid.list[[8]])),7,8)
    sample.dist[8,]<-c(dist(rbind(sample.centroid.list[[7]],sample.centroid.list[[9]])),7,9)
    sample.dist[9,]<-c(dist(rbind(sample.centroid.list[[8]],sample.centroid.list[[9]])),8,9)
    
    #lets put the cluster averages in a set too
    cluster.avg<-array(NA,dim=c(9,1))
    cluster.avg[1]<-mean(md2$mappingData[1:row.pos[1],1:col.pos[1]])
    cluster.avg[4]<-mean(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),1:col.pos[1]])
    cluster.avg[7]<-mean(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),1:col.pos[1]])
    
    cluster.avg[2]<-mean(md2$mappingData[1:row.pos[1],col.pos[1]:(col.pos[1]+col.pos[2])])
    cluster.avg[5]<-mean(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),col.pos[1]:(col.pos[1]+col.pos[2])])
    cluster.avg[8]<-mean(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),col.pos[1]:(col.pos[1]+col.pos[2])])
    
    cluster.avg[3]<-mean(md2$mappingData[1:row.pos[1],(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    cluster.avg[6]<-mean(md2$mappingData[row.pos[1]:(row.pos[1]+row.pos[2]),(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    cluster.avg[9]<-mean(md2$mappingData[(row.pos[1]+row.pos[2]):(row.pos[1]+row.pos[2]+row.pos[3]),(col.pos[1]+col.pos[2]):(col.pos[1]+col.pos[2]+col.pos[3])])
    
    
    ## proportion of DIO in clusters
    # IM.prop<-c(sum(dio.annot[1:col.pos[1],"IM_Cut"])/col.pos[1],
    #            sum(dio.annot[(col.pos[1]+1):(col.pos[1]+col.pos[2]),"IM_Cut"])/col.pos[2],
    #            sum(dio.annot[(col.pos[1]+col.pos[2]+1):(col.pos[1]+col.pos[2]+col.pos[3]),"IM_Cut"])/col.pos[3])
    
    ## tissue partioning
    column_split <- factor(c(rep("(1)", md2$col.size[md2$col.position["left"]]),
                             rep("(2)", md2$col.size[md2$col.position["middle"]]),
                             rep("(3)", md2$col.size[md2$col.position["right"]])), 
                           levels= c("(1)", "(2)", "(3)"))
    tissue_part<-cbind(tissueannot,column_split)
    tissue_part_summary<-table(tissue_part$tissue,tissue_part$column_split)
  }
  cluster_list <- list("cluster_plot" = clusterplot, 
                       #"DIO_prop" = IM.prop, 
                       "gene_dist" = gene.dist, 
                       "sample_dist" = sample.dist, 
                       "row_size"=row.pos, 
                       "col_size"=col.pos, 
                       "cluster_avg"=cluster.avg, 
                       "tissue_part_summary"=tissue_part_summary, 
                       "tissue_annot"=tissueannot,
                       "gene_annot"=geneannot,
                       "clustered_data"=clustered_data)
  return(cluster_list)
}
