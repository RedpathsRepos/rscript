# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Create plot of  
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("/home/miles/RNASeq/")
setwd("/home/miles/Dropbox/Papers/402_RNASeq/working")

#======================================================================================================================#

dat<- read.table("mds.data.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)

#windows(25,25)
par(mfcol=c(1,1),pin=c(5,4),cex.lab=1.1,cex=1.25,lwd=2,oma=c(3,2,0,0))

plot(-22,-22,xlim=c(-1,1),ylim=c(-1 ,1),xlab="Dimension 1", ylab="Dimension 2",pch=24, col="black", bg="white",cex=0.00001,cex.axis=.95)
abline(h=0,lwd=.1)
abline(v=0,lwd=.1)


hatch=dat[36:70,]
points(hatch[,1],hatch[,2],pch=21,bg="green",col="black",cex=2)
WW=dat[1:35,]
points(WW[,1],WW[,2],pch=21,bg="blue",col="black",cex=2)






key=c("Wild fish","F1 hatchery fish")
#legend(-1,1, key, pch=c(19),col=c("blue","green"),cex=0.8)
#legend(-1,1, key, pch=c(21),pt.bg=c("blue","green"),cex=0.65)

#siblings

points(hatch[3,1],hatch[3,2],pch=21,bg="red")
points(hatch[10,1],hatch[10,2],pch=21,bg="red")
points(hatch[16,1],hatch[16,2],pch=21,bg="red")
points(hatch[23,1],hatch[23,2],pch=21,bg="red")
points(hatch[4,1],hatch[4,2],pch=21,bg="orange")
points(hatch[6,1],hatch[6,2],pch=21,bg="orange")
points(hatch[14,1],hatch[14,2],pch=21,bg="orange")
points(hatch[19,1],hatch[19,2],pch=21,bg="orange")
points(hatch[5,1],hatch[5,2],pch=21,bg="blue")
points(hatch[11,1],hatch[11,2],pch=21,bg="blue")
points(hatch[18,1],hatch[18,2],pch=21,bg="blue")
points(hatch[21,1],hatch[21,2],pch=21,bg="blue")
points(hatch[1,1],hatch[1,2],pch=21,bg="yellow")
points(hatch[8,1],hatch[8,2],pch=21,bg="yellow")
points(hatch[22,1],hatch[22,2],pch=21,bg="yellow")
points(hatch[24,1],hatch[24,2],pch=21,bg="yellow")
points(hatch[2,1],hatch[2,2],pch=21,bg="green")
points(hatch[12,1],hatch[12,2],pch=21,bg="green")
points(hatch[15,1],hatch[15,2],pch=21,bg="green")
points(hatch[20,1],hatch[20,2],pch=21,bg="green")
points(hatch[7,1],hatch[7,2],pch=21,bg="gray")
points(hatch[9,1],hatch[9,2],pch=21,bg="gray")
points(hatch[13,1],hatch[13,2],pch=21,bg="gray")
points(hatch[17,1],hatch[17,2],pch=21,bg="gray")

points(dat[37,1],dat[37,2],pch=22,bg="red")
points(dat[44,1],dat[44,2],pch=22,bg="red")
points(dat[30,1],dat[30,2],pch=22,bg="red")
points(dat[31,1],dat[31,2],pch=22,bg="red")

points(dat[43,1],dat[43,2],pch=22,bg="orange")
points(dat[48,1],dat[48,2],pch=22,bg="orange")
points(dat[25,1],dat[25,2],pch=22,bg="orange")
points(dat[34,1],dat[34,2],pch=22,bg="orange")

points(dat[38,1],dat[38,2],pch=22,bg="blue")
points(dat[45,1],dat[45,2],pch=22,bg="blue")
points(dat[26,1],dat[26,2],pch=22,bg="blue")
points(dat[35,1],dat[35,2],pch=22,bg="blue")

points(dat[40,1],dat[40,2],pch=22,bg="yellow")
points(dat[41,1],dat[41,2],pch=22,bg="yellow")
points(dat[29,1],dat[29,2],pch=22,bg="yellow")
points(dat[32,1],dat[32,2],pch=22,bg="yellow")

points(dat[42,1],dat[42,2],pch=22,bg="green")
points(dat[46,1],dat[46,2],pch=22,bg="green")
points(dat[27,1],dat[27,2],pch=22,bg="green")
points(dat[33,1],dat[33,2],pch=22,bg="green")

points(dat[39,1],dat[39,2],pch=22,bg="gray")
points(dat[47,1],dat[47,2],pch=22,bg="gray")
points(dat[28,1],dat[28,2],pch=22,bg="gray")
points(dat[36,1],dat[36,2],pch=22,bg="gray")
