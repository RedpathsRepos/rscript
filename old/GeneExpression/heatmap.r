setwd("C:/POPS")
counts<- read.table("counts.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
counts2=as.matrix(counts[1:60,2:25])

counts3=counts2[,-c(3,7,20)]
counts3=counts3[-c(1:6,8,12,15:16,19,24:27,37:42),]
colnames(counts3)<-c("Wild","Wild","Wild","Wild","Wild","Wild","Wild","Wild","Wild","Wild","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery","Hatchery")
rownames(counts3)<-1:39
heatmap(t(counts3),Rowv=NA,Colv=NA,col = heat.colors(160), scale="column", margins=c(5,10),xlab="Gene Ids",ylab="Fish Type")
heatmap(t(counts3),Rowv=NA,Colv=NA,col = cm.colors(160), scale="column", margins=c(5,10),xlab="Gene Ids",ylab="Fish Type")

heatmap(t(counts2),Rowv=NA,Colv=NA,col = cm.colors(16), scale="column", margins=c(5,10))  #terrain.colors, rainbow, heat.colors, cm.colors, topo.colors