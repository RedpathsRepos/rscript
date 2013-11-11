# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Performs differetial gene expression usinge edgeR 
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("/home/miles/Dropbox/Papers/402_RNASeq/working")
library(DESeq2)
source("http://bioconductor.org/biocLite.R")
biocLite("pasilla")
library("pasilla")
data("pasillaGenes")
countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,levels=c("untreated","treated"))
dds



#======================================================================================================================#
dat.all <- read.csv("genecounts.csv", sep = '\t', header=TRUE, row.names=1)
sam <- cbind(1:length(colnames(dat.all)), colnames(dat.all))      # use this to determine whch fish to use
sam2 <- cbind(sam, substr(sam[,2],16,20), substr(sam[,2],22,22), substr(sam[,2],24,24))
sam2
samples <- sam2[which(sam2[, 4] == "M"), ]           # isolate by sex of offspring
samples <- samples[which(samples[, 5] != "S"), ]           # isolate by matrix
dat <- dat.all[, as.numeric(as.character(samples[, 1]))]                            # index all wxw and hxh females 



dat.names <- colnames(dat)
dat.names
group <- cbind(substr(dat.names, 16, 20),substr(dat.names, 22, 24))
group
colnames(group) <- c("condition","type")
rownames(group) <- substr(dat.names, 1, 8)
colnames(dat) <- substr(dat.names, 1, 8)
dds <- DESeqDataSetFromMatrix(countData = dat,colData = as.data.frame(group),design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition, levels=c("WFxWM", "HFxHM"))
dds$condition
dds

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
head(res)
plotMA(dds)



rld <- rlogTransformation(dds, blind=TRUE)   #takes a long time!  up to an hour or more!
vsd <- varianceStabilizingTransformation(dds)

library("RColorBrewer")
library("gplots")

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, type, sep = " : "))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin = c(13, 13))

print(plotPCA(rld, intgroup = c("condition", "type")))
print(plotPCA(rld, intgroup = c("condition")))
print(plotPCA(rld, intgroup = c("condition"), ntop = 100))

