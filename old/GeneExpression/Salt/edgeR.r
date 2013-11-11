#source("http://bioconductor.org/biocLite.R")       #run these two lines to installl package
#biocLite("edgeR")
library(edgeR)    
setwd("C:/POPS/Salt/Counted_files")
dat<- read.table("GeneCountsPCS_hour24.counted", header=FALSE, sep=" ", na.strings="?",dec=".", strip.white=TRUE,row.names=1)


group=factor(c("control","control","control","salt","salt","salt"))
colnames(dat)<-c("control1","control2","control3","salt1","salt2","salt3")

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
top <- topTags(et,n=2000) #take top n genes
write.table(top,file="PK24hr.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

