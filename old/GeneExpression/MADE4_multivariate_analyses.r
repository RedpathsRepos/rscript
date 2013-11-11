#source("http://bioconductor.org/biocLite.R")  #   install it
#biocLite("made4")
#go here for walkthrough: http://www.bioconductor.org/packages/release/bioc/vignettes/made4/inst/doc/introduction.pdf

setwd("C:/POPS/made4")  #Set the working directory.
counts<- read.table("GeneCounts2.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) #males

library(made4)
library(ade4)

row.names(counts)<-counts[,1]
counts2=array2ade4(counts[,-1],pos=TRUE)
overview(counts2)

k.classes=c("wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch")
counts3=ord(counts2, type="coa", classvec=k.classes)
plot(counts3)
summary(counts3)
plot(counts3,classvec=k.classes,arraycol=c("red","blue"))

plotgenes(counts3)
plotgenes(counts3,nlab=2,col="red")
plotarrays(counts3, arraylabeles=k.classes)

counts4=counts3$ord
plot.ord(counts3, axes1=1, axes2=2, arraycol=NULL, genecol="gray25", nlab=10, genelabels= NULL, classvec=NULL) #works if you get rid of all arguments but the first
#counts3=ord(counts2, type="pca", classvec=NULL)

################################################################################
#Detailed plots
n1=1
n3=2
plot(counts3$ord$co[,n1],counts3$ord$co[,n3],xlab="Axis 1",ylab="Axis 2")

abline(0,0)
lines(c(0,0),c(-10,10))
points(counts3$ord$co[1:12,n1],counts3$ord$co[1:12,n3],pch=21,bg="blue",cex=2)
points(counts3$ord$co[13:24,n1],counts3$ord$co[13:24,n3],pch=21,bg="green",cex=2)
#text(counts3$ord$co[1,n1],counts3$ord$co[1,n3],"4/28",pos=1)

key=c("Hatchery","Wild")
legend(-.37,0.31, key,pch=21,pt.bg=c("blue","green"))
########################################################################
#Between Group Analysis

k.bga<-bga(counts2, tupe="coa", classvec=k.classes)
plot.bga(k.bga)
between.graph(k.bga, ax=1)

##top genes
top=topgenes(counts2)
a=which(counts[,1]==top[1])
counts[a,]