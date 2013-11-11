# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Performs differetial gene expression usinge edgeR 
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("/home/miles/RNASeq/")
setwd("C:/Dropbox/InPrep/RNASeq/working/DifferentialExpression")
library(edgeR)

#======================================================================================================================#
#dat.all <- read.csv("genecounts.consolidated.csv", sep = '\t', header=TRUE, row.names=1)
dat.all <- read.csv("genecounts.csv", sep = '\t', header=TRUE, row.names=1)
sam <- cbind(1:length(colnames(dat.all)), colnames(dat.all))                             # use this to determine whch fish to use
sam2 <- cbind(sam, substr(sam[,2],16,20), substr(sam[,2],22,22), substr(sam[,2],24,24))
sam2
samples <- sam2
samples <- samples[which(samples[, 3] != "HFxWM"), ]                                     # isolate by matrix type
samples <- samples[which(samples[, 3] != "WFxHM"), ]                                     # isolate by matrix type
dat <- dat.all[, as.numeric(as.character(samples[, 1]))]                                 # index all wxw and hxh females 

dat.names <- colnames(dat)
dat.names
group <- substr(dat.names, 16, 20)

y <- DGEList(counts=dat,group=group)
dim(y)                                 # total number of contigs and sample
y$samples                              # total number of reads per sample
levels(y$samples$group)                # check that all groups are represented

#filtering
filter <- rowSums(cpm(y) > 1) >= 35   # keep genes with at least 1 count per million reads in at least 4 samples
y.filtered <- y[filter,]
dim(y.filtered)                       # consider removing the outlier individuals samples 10 and 22

#normalize
y.normalized <- calcNormFactors(y.filtered)

#data exploration
barplot(y.normalized$samples$lib.size*1e-6, ylab = "library size (millions)")   # creates barplot of filtered read #
plotMDS(y.normalized, labels=group)                                             # label options group for group, dat.names for individual sample names
#s=plotMDS(y.normalized)                                                        # get coordinates for MDS

dat.names <- colnames(dat)
dat.names
design <- cbind(substr(dat.names, 16, 20),substr(dat.names, 22, 24))
colnames(design) <- c("treatment","matrix")
rownames(design) <- substr(dat.names, 1, 8)
head(design)
matrix <- factor(design[, 2])
treatment <- factor(design[, 1])

design <- model.matrix(~matrix + treatment)

# estimate the dispersion
y.processed <- estimateGLMCommonDisp(y.normalized, design, verbose = TRUE)   # estimate the overall dispersion
y.processed <- estimateGLMTrendedDisp(y.processed, design)                   # estimate gene-wise dispersion rates
y.processed <- estimateGLMTagwiseDisp(y.processed, design)
plotBCV(y.processed)

# fit the genewise glms
fit <- glmFit(y.processed, design)
lrt <- glmLRT(fit)

topTags(lrt)
summary(de <- decideTestsDGE(lrt))
top.names <- rownames(topTags(lrt, 200))    # select the top xxx genes based on FDR p-value
m1 <- match(top.names, rownames(dat))
dat <- dat[m1,]                             # create data set with only DE genes

top <- topTags(lrt, 200)
write.table(top,file="topgenes_new.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)
asp <- cpm(y.processed)[rownames(top),]

plotSmear(lrt, de.tags = rownames(top))
abline(h=c(-1,1), col = "blue")   #lines equal 2 fold up or down

#recreate mds plot with top genes
y <- DGEList(counts=dat,group=group)
dim(y)                                 # total number of contigs and sample
y$samples                              # total number of reads per sample
levels(y$samples$group)                # check that all groups are represented
filter <- rowSums(cpm(y) > 1) >= 35    # keep genes with at least 1 count per million reads in at least 4 samples
y.filtered <- y[filter,]
dim(y.filtered)                       # consider removing the outlier individuals samples 10 and 22
#normalize
y.normalized <- calcNormFactors(y.filtered)
plotMDS(y.normalized, labels=group, top = 150) 
points.mds=plotMDS(y.normalized, labels = group, top = 150)  
points.data <- cbind(points.mds$x, points.mds$y)
write.table(points.data,file="mds.data.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)


