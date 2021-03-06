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
#samples <- samples[which(samples[, 4] == "F"), ]                     # isolate by sex of offspring
#samples <- samples[which(samples[, 3] != "HFxWM"), ]           # isolate by matrix type
#samples <- samples[which(samples[, 3] != "WFxHM"), ]           # isolate by matrix type
samples <- samples[which(samples[, 5] == "K"), ]               # isolate by specific matrix

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

y <- estimateCommonDisp(y.normalized, verbose = TRUE)
y <- estimateTagwiseDisp(y)
plotBCV(y)

# perform exact test (only valid for 2 treatments) - most basic of tests
basic <- exactTest(y)
basic2 <- basic$table  # note that these have not been corrected for False Discovery
hist(basic2[, 3])  # distribution of all pvalues                                            
summary(de <- decideTestsDGE(basic))

basic2 <- topTags(basic, n = 200)
basic2 <- basic2$table
write.table(basic2, file="TopGenes_ExactTest_Females.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)  # top 500 genes


# begin GLm approach
design <- cbind(substr(dat.names, 16, 20),substr(dat.names, 22, 24))  # using the glm approach, requires a design matrix
colnames(design) <- c("treatment","matrix")
rownames(design) <- substr(dat.names, 1, 8)
head(design)
matrix <- factor(design[, 2])
treatment <- factor(design[, 1])
treatment <- relevel(treatment, ref = "WFxWM")

y <- DGEList(counts=dat,group=treatment)
y.normalized <- calcNormFactors(y.filtered)                   # normalize
plotMDS(y.normalized)

design <- model.matrix(~treatment)
rownames(design) <- colnames(y.normalized)

# estimate the dispersion
y.processed <- estimateGLMCommonDisp(y.normalized, design, verbose = TRUE)  # estimate the overall dispersion
y.processed <- estimateGLMTrendedDisp(y.processed, design)                  # estimate gene-wise dispersion rates
y.processed <- estimateGLMTagwiseDisp(y.processed, design)
plotBCV(y.processed)

# fit the genewise glms
fit <- glmFit(y.processed, design)
lrt <- glmLRT(fit)

summary(de <- decideTestsDGE(lrt))


top <- topTags(lrt, n = 200)
write.table(top,file="Topgenes_GLM_BothSexes.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)


# rerun with blocking
design <- model.matrix(~matrix+treatment)
rownames(design) <- colnames(y.normalized)

# estimate the dispersion
y.processed <- estimateGLMCommonDisp(y.normalized, design, verbose = TRUE)  # estimate the overall dispersion
y.processed <- estimateGLMTrendedDisp(y.processed, design)                  # estimate gene-wise dispersion rates
y.processed <- estimateGLMTagwiseDisp(y.processed, design)
plotBCV(y.processed)

# fit the genewise glms
fit <- glmFit(y.processed, design)
lrt <- glmLRT(fit)

summary(de <- decideTestsDGE(lrt))
top <- topTags(lrt, n = 200)
write.table(top,file="Topgenes_GLMBlocked_Females.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)










#check whether there is a genuine need for blocking, by making appropriate constrasts
#make contrasts from design matrix
head(design)
# IMPORTANT : HxH is "contr.treatment" 
WvsH <- makeContrasts(treatmentHFxHM-treatmentWFxWM, levels = design)
WvsH <- makeContrasts(treatmentWW, levels = design)
hybrids <- makeContrasts(treatmentWH-treatmentHW, levels = design)
WvsH <- makeContrasts(treatmentWFxWM-Intercept, levels = design)

#can also make more interesting contrasts (i.e. 1 contrast vs average of 2 contrasts)


# run likelihood ratio tests with appropriate contrasts
lrt <- glmLRT(fit, contrast = WvsH)  #make pair-wise comparisons
topTags(lrt, n = 100)
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


#========================================================================================================================#
#Extra code:

#s=plotMDS(y.normalized)                                      # get coordinates for MDS
#asp=cbind(s$x,s$y)
#write.table(asp,file="hybrids_MDScoord.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)

#counts <- cpm(y.normalized, normalized.lib.sizes=TRUE)
#pc <- prcomp(counts)
#plot(pc$rotation[1:23, 1], pc$rotation[1:23, 2], pch=21, bg="orange")
#points(pc$rotation[24:45, 1], pc$rotation[24:45, 2], pch=21, bg="blue")




