#source("http://bioconductor.org/biocLite.R")     install it
#biocLite("made4")



setwd("C:/POPS/made4")  #Set the working directory.
counts<- read.table("Male_GeneCounts_v1_f_best_m5_k1_p.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) #males
countsf<- read.table("Female_GeneCounts_v1_f_best_m5_k1_p.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)        #females

counts=cbind(counts,countsf[,-1])
colnames(counts) <- c("Contig","MaleW","MaleW","MaleW","MaleW","MaleH","MaleH","MaleH","MaleH","FemaleH","FemaleH","FemaleH","FemaleH","FemaleH","FemaleW","FemaleW")
colnames(counts) <- c("Contig","W5/17","W5/17","W5/24","W5/24","H5/10","H4/28","H5/17","H2009")
colnames(counts) <- c("Contig","H4/28","H5/24","H5/24","H2009","H2009","W4/28","W4/28")

library(made4)
counts2=array2ade4(counts[,-1])
counts3=ord(counts2, type="coa", classvec=NULL)
plot(counts3)
#plot.ord(counts3, axes1=1, axes2=2, arraycol=NULL, genecol="gray25", nlab=10, genelabels= NULL, classvec=NULL) #works if you get rid of all arguments but the first
#counts3=ord(counts2, type="pca", classvec=NULL)
#################################################################################
#Females
n1=1
n3=2
plot(counts3$ord$co[,n1],counts3$ord$co[,n3],xlab="Axis 1",ylab="Axis 2")

abline(0,0)
lines(c(0,0),c(-10,10))
points(counts3$ord$co[1:5,n1],counts3$ord$co[1:5,n3],pch=21,bg="blue",cex=2)
points(counts3$ord$co[6:7,n1],counts3$ord$co[6:7,n3],pch=21,bg="green",cex=2)
text(counts3$ord$co[1,n1],counts3$ord$co[1,n3],"4/28",pos=1)
text(counts3$ord$co[2,n1],counts3$ord$co[2,n3],"5/24",pos=1,col="red")
text(counts3$ord$co[3,n1],counts3$ord$co[3,n3],"5/24",pos=1,col="red")
text(counts3$ord$co[4,n1],counts3$ord$co[4,n3],"2009",pos=2)
text(counts3$ord$co[5,n1],counts3$ord$co[5,n3],"2009",pos=1)
text(counts3$ord$co[6,n1],counts3$ord$co[6,n3],"4/28",pos=3,col="purple")
text(counts3$ord$co[7,n1],counts3$ord$co[7,n3],"4/28",pos=1,col="purple")
key=c("Hatchery","Wild")
legend(-.37,0.31, key,pch=21,pt.bg=c("blue","green"))

########################################################################
#Males
 
n1=1
n3=2
plot(counts3$ord$co[,n1],counts3$ord$co[,n3],xlab="Axis 1",ylab="Axis 2")

abline(0,0)
lines(c(0,0),c(-10,10))
points(counts3$ord$co[1:4,n1],counts3$ord$co[1:4,n3],pch=21,bg="blue",cex=2)
points(counts3$ord$co[5:8,n1],counts3$ord$co[5:8,n3],pch=21,bg="green",cex=2)
text(counts3$ord$co[1,n1],counts3$ord$co[1,n3],"5/17",pos=3)
text(counts3$ord$co[2,n1],counts3$ord$co[2,n3],"5/17",pos=1)
text(counts3$ord$co[3,n1],counts3$ord$co[3,n3],"5/24",pos=3)
text(counts3$ord$co[4,n1],counts3$ord$co[4,n3],"5/24",pos=2)
text(counts3$ord$co[5,n1],counts3$ord$co[5,n3],"5/10",pos=1)
text(counts3$ord$co[6,n1],counts3$ord$co[6,n3],"4/28",pos=4)
text(counts3$ord$co[7,n1],counts3$ord$co[7,n3],"5/17",pos=1)
text(counts3$ord$co[8,n1],counts3$ord$co[8,n3],"2009",pos=1)
key=c("Hatchery","Wild")
legend(.31,0.37, key,pch=21,pt.bg=c("blue","green"))

colnames(counts) <- c("Contig","W5/17","W5/17","W5/24","W4/20","H5/10","H4/28","H5/17","H2009")