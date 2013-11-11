#This script takes data from one catogory of multiple rows and places all like categories into multiple columns

setwd("C:/Users/christim/Desktop/Cluster")
All<- read.table("AllOut.txt", header=F, sep="\t", na.strings="?", dec=".", strip.white=TRUE)

Ids=unique(All[,2])


func=function(one)
{
Id2=Ids[z]
gos=All[which(All[,2]==Id2),]
out1=t(gos)
out2=c(out1[3,],out1[4,])
out2=t(out2)
out3=cbind(as.character(Id2),out2)

write.table(out3,file="AllGOs.txt",col.names=F,sep="\t",append=T)
}
C1<- for(z in (1:length(Ids))) sapply(z,func)

