
#blocking by matrix

library(edgeR)

setwd("~/R/POPS/GeneExpression")
dat<- read.table("GeneCounts_males.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE,row.names=1)
group=factor(c("hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild"))


#using getPriorN as recommended in vijay Molec Ecol paper
y <- DGEList(counts=dat,group=group)
colnames(y)
dim(y)        #total number of unique tags
y$samples     #number of tags ("genes") per sample


#create targets for blocking

targets1 <- colnames(dat)
targets2 <- c("W","J","T","I","W","L","J","T","L","S","I","S","I","L","W","S","S","T","J","I","L","T","W","J")
targets3 <- group

targets <- as.data.frame(cbind(targets1,targets2,as.character(targets3)))
colnames(targets) <- c("FileName", "Litter", "Treatment")
targets <- targets[order(targets[,2], targets[,3]),]             #sort by columns 1 and 2  ; may need to omit this step (file order altered)

Litter <- factor(targets$Litter)
Treatment <- factor(targets$Treatment)
design <- model.matrix(~Litter+Treatment)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=7)

topTags(lrt)

summary(lrt <- decideTestsDGE(lrt))    #-1 equals upregulate, 1 equals down regulated:  need to be carefuul though as order (control,vs treatment matters)
top <- topTags(et,n=40) #take top n genes
write.table(top,file="DEgenesboth2.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")






