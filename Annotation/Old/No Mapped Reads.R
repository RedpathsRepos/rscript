# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Uses gene counts to remove find uneccessary contigs because the number of mapped reads are so low
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