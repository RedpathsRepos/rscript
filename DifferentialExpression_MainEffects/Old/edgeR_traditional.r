#==================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script:  Uses edgeR to test for differential gene expression
# Usage notes:  NA
#==================================================================================================================#
# Source files, import packages, set working directory, initialize variables
# source("http://bioconductor.org/biocLite.R")                                  # run these two lines to install package from bioconductor, if not installed
# biocLite("edgeR")
library(edgeR)
directory <- "C:/Dropbox/Papers/402_RNASeq/data" 
input <- "GeneCounts_males2.txt"
output <- "DE_genes.txt"


#==================================================================================================================#
setwd(directory)
dat <- read.table(input, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE, row.names=1)

group <- factor(c("hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch", "hatch","hatch","hatch","hatch",
                  "wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild"))               # convert columns to factors



#using getPriorN as recommended in vijay Molec Ecol paper
y <- DGEList(counts=dat,group=group)
dim(y)

#filter
keep <- rowSums(cpm(y)) >= 2
y <- y[keep,]
dim(y)

#recompute library sizes
y$samples$lib.size <- colSums(y$counts)

#normalize
y <- calcNormFactors(y)
y$samples

#plot MDS of treatments
plotMDS(y)

#s=plotMDS(y)                    #get coordinates for MDS
#asp=cbind(s$x,s$y)
#write.table(asp,file="Males_MDScoord.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)


#estimate dispersion
#prior=getPriorN(y)     #reocmmended in Molec Ecol paper, not sure how to implement

 y <- estimateCommonDisp(y, verbose=TRUE)
 y <- estimateTagwiseDisp(y)
# y <- estimateTagwiseDisp(y,prior.df=prior)    #I am not sure if this is the correct implementaion
 plotBCV(y)
 
#test for DE
et <- exactTest(y)
top <- topTags(et)
top
cpm(y)[rownames(top), ]

#total number of genes at 5% FDR
summary(de <- decideTestsDGE(et))    #-1 equals upregulate, 1 equals down regulated:  need to be carefuul though as order (control,vs treatment matters)
top <- topTags(et,n=120) #take top n genes
write.table(top,file="DEgenesMales.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

