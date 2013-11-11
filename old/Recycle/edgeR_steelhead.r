#source("http://bioconductor.org/biocLite.R")       #run these two lines to installl package
#biocLite("edgeR")
library(edgeR)

setwd("~/R/POPS/GeneExpression")
dat<- read.table("GeneCounts_males.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE,row.names=1)
group=factor(c("hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild"))

#Two quick examples

y <- DGEList(counts=dat,group=group)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)

write.table(et$table,file="Qvalues.txt",col.names=TRUE,sep="\t",append=FALSE)


design <- model.matrix(~group)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt,n=200)


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

