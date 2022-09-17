#the code below calls this function a few times

convert_to_numeric<-function(convertdata=NULL) {
  
  convertdata<-as.data.frame(convertdata)
  
  therownames<-rownames(convertdata)
  
  convertdata <- apply(convertdata, 2,function(x) as.numeric(as.character(x)))
  
  rownames(convertdata)<-therownames
  
  return(convertdata)
  
}



#load the drug sensitivity data

drugScreen <- read.csv("primary-screen-replicate-collapsed-logfold-change (1).csv")

colnames(drugScreen) <- substring(colnames(drugScreen), 5, 13)

drugScreen <- drugScreen[-which(as.numeric(regexpr("FAILED", drugScreen[,1]) )>0),]

rownames(drugScreen) <- drugScreen[,1]

drugScreen <- drugScreen[,-1]

dim(drugScreen) #568 4686



## limit to training cell set

cellannot<-read.csv("ccle_annotation 12-6-21.csv", header=FALSE,na.strings = "NA",stringsAsFactors = FALSE)

cellannot<-as.data.frame(cellannot[-1,])

colnames(cellannot)<-cellannot[1,]

cellannot<-cellannot[-1,]

drugScreen<-drugScreen[rownames(drugScreen) %in% cellannot[cellannot$test_set_balanced==0,"sampleid"],]

dim(drugScreen) #371 4686



## load the cell line data

cellexp<-readRDS("CCLE_expression.rds")

cellexp<-cellexp[rownames(drugScreen),] #limit to the training cell lines also in the drug data

dim(cellexp) #371 19177



## correlate drugs with tnbctype for the cell lines then choose drugs with any corr < -0.1 (remember we want the negative, inhibitory values)

CellLineByClass<-cellannot[,c("IM_Corr","MSL_Corr","M_Corr")]

rownames(CellLineByClass)<-cellannot$sampleid

CellLineByClass<-CellLineByClass[rownames(drugScreen),] #limit to shared training samples



CellLineByClass<-CellLineByClass[order(rownames(CellLineByClass)),]  #lets order things to be safe

drugScreen<-drugScreen[order(rownames(drugScreen)),]  #lets order things to be safe

CellLineByClass<-convert_to_numeric(CellLineByClass)

drugScreen<-convert_to_numeric(drugScreen)



drugtype <- data.frame(IM_Corr=double(ncol(drugScreen)),
                       
                       MSL_Corr=double(ncol(drugScreen)),
                       
                       M_Corr=double(ncol(drugScreen)),
                       
                       stringsAsFactors=FALSE)

rownames(drugtype)<-colnames(drugScreen)

for (n in 1:nrow(drugtype)) {
  
  drugtype[n,1] <- cor(drugScreen[,n], CellLineByClass[,1],method = "pearson",  use = "complete.obs")
  
  drugtype[n,2] <- cor(drugScreen[,n], CellLineByClass[,2],method = "pearson",  use = "complete.obs")
  
  drugtype[n,3] <- cor(drugScreen[,n], CellLineByClass[,3],method = "pearson",  use = "complete.obs")
  
}



strongdrugs<-data.frame(drugtype[rowSums(drugtype < -0.1) >= 1,])

dim(strongdrugs) #864 3



#grab our selected gene set

canonical_genes<-read.csv("combined top canon sets v5.csv", header=TRUE,na.strings = "NA",stringsAsFactors = FALSE)

canonical_genes<-canonical_genes[,1]

## lets add TUSC2

canonical_genes<-c(canonical_genes,"TUSC2")



## create data set of canon genes and strong drugs on training samples

drugscreenlimit<-drugScreen[,rownames(strongdrugs)]

sharedcells<-intersect(rownames(drugscreenlimit),rownames(cellexp))

length(sharedcells) #371

cellexp<-cellexp[sharedcells,which(colnames(cellexp) %in% canonical_genes)] #limit to the training cell lines also in the drug data, and the selected gene list

drugscreenlimit<-drugscreenlimit[sharedcells,]

dim(drugscreenlimit) #371, 864

dim(cellexp) #371, 2460



cellexp<-cellexp[order(rownames(cellexp)),]  #lets order things to be safe

drugscreenlimit<-drugscreenlimit[order(rownames(drugscreenlimit)),]  #lets order things to be safe



library(ggplot2)

library(matrixStats)

ggplot() + aes(x=(rowSds(t(cellexp)))) +
  
  geom_density() + ggtitle("cell line gene expression") +
  
  xlab("gene sds") +xlim(c(-1,3))

## a value of 0.3 still looks like to cut off a low peak

tcellexp<-t(cellexp)

cellexp_sel<-t(tcellexp[rowSds(tcellexp)>0.3,])

dim(cellexp_sel) ##371 2196



ggplot() + aes(x=(rowSds(t(drugscreenlimit),na.rm = TRUE))) +
  
  geom_density() + ggtitle("compound sensitivity") +
  
  xlab("drug sds") +xlim(c(-0,2))

## that still looks OK



## scale each gene/drug

scellexp<-scale(cellexp_sel)



## scale each sample

scellexp<-t(scale(t(scellexp)))



## since a negative compound sensitivity number means that the growth of the

# cell line is inhibited by the compound, for the heatmap I should plot the inverse

# of that.  Otherwise the genes that map near a drug are ones that high expression correlates

# with decreased sensitivity.

drugscreenlimitneg<-drugscreenlimit*(-1)

sdrugscreenlimitneg<-scale(drugscreenlimitneg)



## scale each sample

sdrugscreenlimitneg<-t(scale(t(sdrugscreenlimitneg)))



dgmapDataneg<-cbind(scellexp,sdrugscreenlimitneg)

dgmapDataneg<-t(dgmapDataneg)

dim(dgmapDataneg) #3060, 371



tissueannotation<-read.csv("ccle_annotation 12-6-21.csv", header=FALSE,na.strings = "NA",stringsAsFactors = FALSE)

genedrugannot<-read.csv("combined top canon sets v5 annotation.csv", header=FALSE,na.strings = "NA",stringsAsFactors = FALSE)

tissueannotation[,1] <- gsub("-", "\\.", tissueannotation[,1])

clustest<-cluster_plus_ultra(mapData=dgmapDataneg,
                             
                             tissueannot=tissueannotation,
                             
                             geneannot=genedrugannot,
                             
                             scale.genes=FALSE, scale.samples=TRUE,
                             
                             annotate_genes=TRUE,
                             
                             cluster.method="sumsquares", col_cluster_method="kmeans", ordering_columns="diagonal",
                             
                             return.plot=TRUE, return.stats=FALSE, return_data=TRUE)



#save the plot as a png file
clustest[["cluster_average"]]
png("heatmapping/cells drugs training 101 adding in TUSC2 neg drugs 3-3-22.png",  3500, 2150)

draw(clustest$cluster_plot,annotation_legend_side = 'left')

dev.off()


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

#save the clustered data and the clustered gene annotation

write.csv(clustest$clustered_data,"various/cells drugs training plus TUCS2 neg drugs 3-2.csv")

write.csv(clustest$gene_annot,"various/cells drugs training plus TUCS2 gene annotation neg drugs 3-2.csv")

install_github("jokergoo/ComplexHeatmap")

write.csv(dgmapDataneg, "dgmapDataneg.csv")
