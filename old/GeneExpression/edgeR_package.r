#source("http://bioconductor.org/biocLite.R")       #run these two lines to installl package
#biocLite("edgeR")
library(edgeR)

setwd("C:/POPS/made4")
dat<- read.table("GeneCounts2.txt", header=FALSE, sep="\t", na.strings="?",dec=".", strip.white=TRUE,row.names=1)
group=factor(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2))
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
