# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Performs differetial gene expression using the edgeR package
# Usage notes:  First uses exact test approach, then uses glm approach
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("/home/miles/RNASeq/")
setwd("C:/Dropbox/InPrep/RNASeq/working/DifferentialExpression")
library(edgeR)


#======================================================================================================================#
dat.all <- read.csv("genecounts.consolidated.csv", sep = '\t', header=TRUE, row.names=1)
sam <- cbind(1:length(colnames(dat.all)), colnames(dat.all))      # use this to determine whch fish to use
samples <- cbind(sam, substr(sam[,2],16,20), substr(sam[,2],22,22), substr(sam[,2],24,24))
samples
#samples <- samples[which(samples[, 4] == "F"), ]              # isolate by sex of offspring
#samples <- samples[which(samples[, 3] != "HFxWM"), ]          # isolate by matrix type
#samples <- samples[which(samples[, 3] != "WFxHM"), ]          # isolate by matrix type
samples <- samples[which(samples[, 5] == "M"), ]               # isolate by specific matrix

dat <- dat.all[, as.numeric(as.character(samples[, 1]))]       # index desired individuals 

dat.names <- colnames(dat)                                     # check to make sure desired individuals were obtained 
dat.names
group <- substr(dat.names, 16, 20)
table(group)

y <- DGEList(counts=dat,group=group)                          # Begin formatting for edgeR
dim(y)                                                        # total number of contigs and sample
y$samples                                                     # total number of reads per sample
levels(y$samples$group)                                       # check that all groups are represented

filter <- rowSums(cpm(y)> .5) >= 2                            # keep genes with at least x count per million reads in at least y samples


y.filtered <- y[filter,]
dim(y.filtered)

y.normalized <- calcNormFactors(y.filtered)                   # normalize

barplot(y.normalized$samples$lib.size*1e-6, ylab = "library size (millions)")  #creates barplot of filtered read #
plotMDS(y.normalized, labels=group)
plotMDS(y.normalized, labels=dat.names, top = 100)            # MDS plot: label options group for group, dat.names for individual sample names

# begin GLm approach
design <- model.matrix(~0+group, data = y.filtered$samples)
colnames(design) <- levels(y.filtered$samples$group)

# estimate the dispersion
y.processed <- estimateGLMCommonDisp(y.normalized, design, verbose = TRUE)  # estimate the overall dispersion
y.processed <- estimateGLMTrendedDisp(y.processed, design)                  # estimate gene-wise dispersion rates
y.processed <- estimateGLMTagwiseDisp(y.processed, design)
plotBCV(y.processed)

fit <- glmFit(y.processed, design)
lrt1 <- glmLRT(fit, contrast=c(0,1,0,-1))   # WW vs HW
lrt2 <- glmLRT(fit, contrast=c(0,0,1,-1))   # WW vs WH

lrt3 <- glmLRT(fit, contrast=c(-1,1,0,0))   # HH vs HW
lrt4 <- glmLRT(fit, contrast=c(-1,0,1,0))   # HH vs WH

summary(de <- decideTestsDGE(lrt1))
summary(de <- decideTestsDGE(lrt2))
summary(de <- decideTestsDGE(lrt3))
summary(de <- decideTestsDGE(lrt4))




# write.table(top,file="Topgenes_GLMBlocked_Females.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)




