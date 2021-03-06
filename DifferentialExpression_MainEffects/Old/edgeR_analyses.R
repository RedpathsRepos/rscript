# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Performs differetial gene expression usinge edgeR 
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("/home/miles/RNASeq/")
setwd("/home/miles/Dropbox/Papers/402_RNASeq/working")
library(edgeR)


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
group <- substr(dat.names, 16, 20)


# Begin formatting for edgeR
y <- DGEList(counts=dat,group=group)
dim(y)  # total number of contigs and sample
y$samples  # total number of reads per sample
levels(y$samples$group)   # check that all groups are represented

#filtering
filter <- rowSums(cpm(y)>1) >= 10 #kepp genes with at least 1 count per million reads in at least 4 samples
y.filtered <- y[filter,]
dim(y.filtered)

#normalize
y.normalized <- calcNormFactors(y.filtered)

#data exploration
barplot(y.normalized$samples$lib.size*1e-6, ylab = "library size (millions)")  #creates barplot of filtered read #
plotMDS(y.normalized, labels=group, top = 100)
#s=plotMDS(y.normalized)                    #get coordinates for MDS
#asp=cbind(s$x,s$y)
#write.table(asp,file="hybrids_MDScoord.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)

#counts <- cpm(y.normalized, normalized.lib.sizes=TRUE)
#pc <- prcomp(counts)
#plot(pc$rotation[1:23, 1], pc$rotation[1:23, 2], pch=21, bg="orange")
#points(pc$rotation[24:45, 1], pc$rotation[24:45, 2], pch=21, bg="blue")

# using the glm approach, requires a design matrix
# read in and format the design matrix
design <- read.table("HWDesignMatrix.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
matrix <- factor(design[, 2])
treatment <- factor(design[, 3])

dat.names <- colnames(dat)
dat.names
group <- cbind(substr(dat.names, 16, 20),substr(dat.names, 22, 24))
group
colnames(group) <- c("condition","type")
rownames(group) <- substr(dat.names, 1, 8)
colnames(dat) <- substr(dat.names, 1, 8)


design <- model.matrix(~0+as.factor(group[, 1]) + as.factor(group[, 2]))
rownames(design) <- colnames(y.normalized)

# estimate the dispersion
y.processed <- estimateGLMCommonDisp(y.normalized, design, verbose = TRUE)  # estimate the overall dispersion
y.processed <- estimateGLMTrendedDisp(y.processed, design)    # estimate gene-wise dispersion rates
y.processed <- estimateGLMTagwiseDisp(y.processed, design)
plotBCV(y.processed)

# fit the genewise glms
fit <- glmFit(y.processed, design)
lrt <- glmLRT(fit)

topTags(lrt)
summary(de <- decideTestsDGE(lrt))
top <- topTags(lrt, 2500)
write.table(top,file="topgenes_males.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)
asp <- cpm(y.processed)[rownames(top),]

#check whether there is a genuine need for blocking, by making appropriate constrasts
#make contrasts from design matrix
head(design)
# IMPORTANT : HxH is "contr.treatment" 
LvsM <- makeContrasts(matrixL-matrixM, levels = design)
WvsH <- makeContrasts(treatmentWW, levels = design)
hybrids <- makeContrasts(treatmentWH-treatmentHW, levels = design)
#can also make more interesting contrasts (i.e. 1 contrast vs average of 2 contrasts)


# run likelihood ratio tests with appropriate contrasts
lrt <- glmLRT(fit, contrast = WvsH)  #make pair-wise comparisons
topTags(lrt)
FDR <- p.adjust(lrt$table$PValue, method = "BH")
sum(FDR < 0.05)

top <- rownames(topTags(lrt, 380))
counts <- cpm(y.processed)[top,]
summary(dt <- decideTestsDGE(lrt))

isDE <- as.logical(dt)
de.names <- rownames(y.processed)[isDE]
de.names  # names of the DE genes

plotSmear(lrt, de.tags = rownames(top))
abline(h=c(-1,1), col = "blue")   #lines equal 2 fold up or down
