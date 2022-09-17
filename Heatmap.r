dgmapDataneg <- read.csv("dgmapDataneg.csv", header=T,row.names = 1, na.strings = "NA", stringsAsFactors = F)
tissueannotation <- read.csv("ccle_annotation 12-6-21.csv", header=F, na.strings = "NA", stringsAsFactors = F)
genedrugannot <- read.csv("combined top canon sets v5 annotation.csv", header = F, na.strings = "NA", stringsAsFactors = F)
#adding more types of annotations
genedrugannot$V4[1] <- 0
genedrugannot$V6[1] <- 1
genedrugannot$V7[1] <- 1
genedrugannot$V8[1] <- 1

annotations <- data.frame(clustest[["gene_annot"]])
tissueannotation[,1]<-gsub("-",'\\.',tissueannotation[,1])

annotations2 <- (annotations$panther_protein[1209:1272])

annotations3 <- data.frame(annotations$panther_protein[2505:2576])
unique(annotations3)
annotations2[1]
for (protein in unique(annotations2)){
  count = length(annotations2[annotations2$V1 == protein])
  print(protein + " : " + count)
}
 clustest<-cluster_plus_ultra_old(mapData=dgmapDataneg,
                                 
                                 tissueannot=tissueannotation,
                                 
                                 geneannot=genedrugannot,
                                 
                                 scale.genes=FALSE, scale.samples=TRUE,
                                 
                                 annotate_genes=TRUE,
                                 
                                 cluster.method="sumsquares", col_cluster_method="kmeans", ordering_columns="diagonal",
                                 
                                 return.plot=TRUE, return.stats=FALSE, return_data=TRUE)

png("large_map.png",  7000, 4300)
 
 
draw(clustest$cluster_plot,annotation_legend_side = 'left')

dev.off()

test123<-tissueannotation[-1,]
colnames(test123)<-tissueannotation[1,]
test123<-test123[-1,]
test123<-test123[,colnames(test123) %in% c("IM_Corr","MSL_Corr","M_Corr")]
test123<-test123[rownames(test123) %in% colnames(dgmapDataneg),]
test123<-test123[order(rownames(test123)),]
therownames<-rownames(test123)
test123$IM_Corr <- lapply(test123$IM_Corr, remove_first_decimal)
test123 <- apply(test123, 2,function(x) as.numeric(as.character(x)))
rownames(test123)<-therownames

test123

remove_first_decimal <- function(string) {
  if (substr(string, 1, 1) == "."){
    string = substr(string, 2, str_length(string))
  }
  return(string)
}

remove_first_decimal('.3')

test1 <- read.csv("test123.csv")
test2 <- read.csv("test1234.csv")

mapdata <- data.frame(clustest[["clustered_data"]]>3)
max(mapdata)



new_df <- data.frame(clustest$clustered_data)

new_df$drug_or_gene <- classifications

new_df = drop_na(new_df)

new_df <- new_df[,1:371]

new_df <- t(new_df)

classification <- annotations$drug_or_gene
#new_df <- new_df[,c(2470:2576)]
#new_df <- new_df[,c(730:808)]

correlations <- cor(new_df, use = "complete.obs")

correlations <- data.frame(correlations)

count = 0
for (i in 1:3060) {
  for (j in 1:3060) {
    if (correlations[i,j] > 0.2 & (classification[i] == "drug" & classification[j] == "gene" | classification[i] == "gene" & classification[j] == "drug") {
      count = count + 1
    }
  }
}

colors <- colorRampPalette(c("blue","white","red"))(99)

png("correlations.png",  7500, 4500)

heatmap(as.matrix(correlations), Colv = NA, Rowv = NA, col = colors, scale="none")

dev.off()
