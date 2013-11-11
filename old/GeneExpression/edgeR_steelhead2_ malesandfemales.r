#source("http://bioconductor.org/biocLite.R")       #run these two lines to installl package
#biocLite("edgeR")
library(edgeR)

setwd("C:/POPS/GeneExpression")
dat<- read.table("GeneCounts_both.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE,row.names=1)
group=factor(c("hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild"))




#using getPriorN as recommended in vijay Molec Ecol paper
y <- DGEList(counts=dat,group=group)
colnames(y)
dim(y)        #total number of unique tags
y$samples     #number of tags ("genes") per sample

#is the number of unique tags greater than the total number of reads in each library (most likely), if so - do some filtering
#filter

keep <- rowSums(cpm(y)>0.5) >= 6         #first number is minimum number of counts per million, if 6 million reads then minimum count will be 6 ; second number is across number of samples < -usually min group size
y <- y[keep,]
dim(y)

#recompute library sizes
y$samples$lib.size <- colSums(y$counts)
y$samples

#normaliez
y <- calcNormFactors(y)
y$samples

#heatmap
#y2 <- predFC(y, prior.count.total=2*ncol(y))  #output can be used for heatmap

#plot MDS of treatments
plotMDS(y)
plotMDS(y, top=300)

#s=plotMDS(y)                    #get coordinates for MDS
#asp=cbind(s$x,s$y)
#write.table(asp,file="Males_MDScoord.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)


#estimate dispersion
#prior=getPriorN(y)     #reocmmended in Molec Ecol paper, not sure how to implement

y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)                      #uses empirical Bayes method, if prior null then df=20
#y <- estimateTagwiseDisp(y,prior.df=5)          #I am not sure if this is the correct implementaion, put more weight on tagwise vs. common
y
plotBCV(y)

#test for DE
et <- exactTest(y)
top <- topTags(et)
top
cpm(y)[rownames(top), ]

#total number of genes at 5% FDR
summary(de <- decideTestsDGE(et))    #-1 equals upregulate, 1 equals down regulated:  need to be carefuul though as order (control,vs treatment matters)
top <- topTags(et,n=1000) #take top n genes
write.table(top,file="DEgenesboth.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
