
setwd("C:/POPS/GeneExpression")
dat<- read.table("Males_MDScoord.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE,row.names=1)



#windows(25,25)
par(mfcol=c(1,1),pin=c(5,4),cex.lab=1.1,cex=1.25,lwd=2,oma=c(3,2,0,0))

plot(-2,-2,xlim=c(-1,1),ylim=c(-1 ,1),xlab="Dimension 1", ylab="Dimension 2",pch=24, col="black", bg="white",cex=0.00001,cex.axis=.95)
abline(h=0,lwd=.1)
abline(v=0,lwd=.1)



hatch=dat[1:12,]
points(hatch[,1],hatch[,2],pch=21,bg="blue",col="black",cex=0.00000000000000000000001)
wild=dat[13:24,]
points(wild[,1],wild[,2],pch=21,bg="green",col="black",cex=0.000000000000000000001)

key=c("Wild fish","F1 hatchery fish")
#legend(-1,1, key, pch=c(19),col=c("blue","green"),cex=0.8)
#legend(-1,1, key, pch=c(21),pt.bg=c("blue","green"),cex=0.65)

#siblings

points(hatch[4,1],hatch[4,2],pch=21,bg="red")
points(hatch[6,1],hatch[6,2],pch=21,bg="yellow")
points(hatch[1,1],hatch[1,2],pch=21,bg="blue")
points(hatch[10,1],hatch[10,2],pch=21,bg="green")
points(hatch[2,1],hatch[2,2],pch=21,bg="orange")
points(hatch[3,1],hatch[3,2],pch=21,bg="gray")
points(hatch[11,1],hatch[11,2],pch=21,bg="red")
points(hatch[9,1],hatch[9,2],pch=21,bg="yellow")
points(hatch[5,1],hatch[5,2],pch=21,bg="blue")
points(hatch[12,1],hatch[12,2],pch=21,bg="green")
points(hatch[7,1],hatch[7,2],pch=21,bg="orange")
points(hatch[8,1],hatch[8,2],pch=21,bg="gray")

hatch=dat
points(hatch[13,1],hatch[13,2],pch=22,bg="red")
points(hatch[14,1],hatch[14,2],pch=22,bg="yellow")
points(hatch[23,1],hatch[23,2],pch=22,bg="blue")
points(hatch[16,1],hatch[16,2],pch=22,bg="green")
points(hatch[19,1],hatch[19,2],pch=22,bg="orange")
points(hatch[18,1],hatch[18,2],pch=22,bg="gray")
points(hatch[20,1],hatch[20,2],pch=22,bg="red")
points(hatch[21,1],hatch[21,2],pch=22,bg="yellow")
points(hatch[15,1],hatch[15,2],pch=22,bg="blue")
points(hatch[17,1],hatch[17,2],pch=22,bg="green")
points(hatch[24,1],hatch[24,2],pch=22,bg="orange")
points(hatch[22,1],hatch[22,2],pch=22,bg="gray")



