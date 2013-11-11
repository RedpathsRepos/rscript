setwd("C:/POPS/")
DE<- read.table("DEgenesMales.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)
GO<- read.table("GOoutputMales.txt", header=FALSE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)

contigs=GO[,1]
OUT=NULL
for (i in 1:length(contigs)) {
  gterm=GO[i,1]
  m1=match(gterm,DE[,1])
  out=cbind(DE[m1,],GO[i,])
  OUT=rbind(OUT,out)}
  
colnames(OUT)<-c("contig","logFC","logCPM","pvalue","FDR","contig.1","GOterm","GOdescription")
write.table(OUT,file="GO_and_pvalues.txt",col.names=TRUE,row.names=FALSE,sep="\t",append=FALSE)