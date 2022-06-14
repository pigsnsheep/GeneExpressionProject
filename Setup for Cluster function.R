#the code below calls this function a few times

convert_to_numeric<-function(convertdata=NULL) {
  
  convertdata<-as.data.frame(convertdata)
  
  therownames<-rownames(convertdata)
  
  convertdata <- apply(convertdata, 2,function(x) as.numeric(as.character(x)))
  
  rownames(convertdata)<-therownames
  
  return(convertdata)
  
}
rownames(drugdata_csv) <- drugdata_csv[,1]
colnames(drugdata_csv) <- substring(colnames(drugdata_csv),5,13)
drugdata_csv <- drugdata_csv[,-c(1)]

drugdata <- drugdata_csv[which(as.numeric(regexpr("FAILED", rownames(drugdata_csv)) == 0),]



