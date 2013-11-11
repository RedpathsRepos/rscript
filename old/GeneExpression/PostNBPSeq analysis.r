setwd("C:/Users/christim/Desktop/Cluster")
qvals<- read.table("Qvalues-conservative.txt", header=F, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
salmongo<- read.table("salmongos.txt", header=F, sep="\t", na.strings="?", dec=".", strip.white=TRUE)

qvals2=as.character(qvals[,1])
salmongo2=as.character(salmongo[,1])

OutH=NULL
func=function(one)
{
contig=qvals2[z]
gos=salmongo[which(salmongo2[]==contig),]
#gos2=if (length(gos[,1])>0) cbind(contig,gos)
OutH<<-rbind(OutH,gos)

}
C1<- for(z in (1:length(qvals[,1]))) sapply(z,func)


General=as.data.frame(table(OutH[,4]))
bb=as.data.frame(table(OutH[,3]))
bc=bb[which(bb[,2]>0),]
Specific=bc[order(bc[,2]),]

write.table(General,file="GoOutGeneral.txt",col.names=T,sep="\t",append=T)
write.table(Specific,file="GoOutSpecific.txt",col.names=T,sep="\t",append=T)


bb1=as.data.frame(table(OutH[,1]))
bc1=bb1[which(bb1[,2]>0),]
Counts=bc1[order(bc1[,2]),]
bd=match(Counts[,1],qvals[,1])
Counts2=cbind(Counts,qvals[bd,])
write.table(Counts2,file="GoOutCounts.txt",col.names=T,sep="\t",append=T)
write.table(OutH,file="AllOut.txt",col.names=T,sep="\t",append=T)
