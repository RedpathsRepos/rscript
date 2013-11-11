library(gplots)
counts <- dat[1:50,]
cor.counts <- cor(counts)

#heatmap.2(cor(counts), Rowv="FALSE")

hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
heatmap.2(cor(counts), col=rev(hmcol),  dendrogram = "none", trace="none", margin = c(13, 13))

counts2 <- counts[,c(1:21,31:42)]
heatmap.2(cor(counts2), col=hmcol, Rowv="FALSE", Colv="FALSE", trace = "none")

#heatmap(t(counts3),Rowv=NA,Colv=NA,col = heat.colors(16), scale="column", margins=c(5,10),xlab="Gene Ids",ylab="Fish Type")
#heatmap.2(t(counts3[1:10,]),Rowv=NA,Colv=NA,col = heat.colors(10), scale="column", margins=c(5,10),xlab="Gene Ids",ylab="Fish Type")
#heatmap(t(counts3[1:20, ]),Rowv=NA,Colv=NA,col = cm.colors(10), scale="column", margins=c(5,10),xlab="Gene Ids",ylab="Fish Type")
#heatmap(t(counts3[1:10, ]),Rowv=NA,Colv=NA,col = cm.colors(16), scale="column", margins=c(5,10))  #terrain.colors, rainbow, heat.colors, cm.colors, topo.colors

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


library("RColorBrewer")
library("gplots")


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
heatmap.2(as.matrix(mat[1:20,]), dendrogram = "none", trace="none", col = rev(hmcol), margin = c(13, 13))