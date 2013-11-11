# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Plots the top genes using MDS for the top X genes using edgeR MDS
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables

setwd("C:/Dropbox/InPrep/RNASeq/working/DifferentialExpression_MainEffects")
library(edgeR)
genes <- read.table("TopGenes_consolidated_mosthits.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)

#======================================================================================================================#

# Make sure the below is identical to samples selected in 'edgeR_MainEffects_PairedTests'
dat.all <- read.csv("genecounts.consolidated.mosthits.csv", sep = '\t', header=TRUE, row.names=1)
sam <- cbind(1:length(colnames(dat.all)), colnames(dat.all))  # use this to determine whch fish to use
samples <- cbind(sam, substr(sam[,2],16,20), substr(sam[,2],22,22), substr(sam[,2],24,24))
samples
samples <- samples[which(samples[, 3] != "HFxWM"), ]          # isolate by matrix type
samples <- samples[which(samples[, 3] != "WFxHM"), ]          # isolate by matrix type
samples <- samples[-c(31, 2, 18),]                            # these samples are extreme outliers via MDS analysis

dat <- dat.all[, as.numeric(as.character(samples[, 1]))]      # index desired individuals 

m1=match(rownames(genes), rownames(dat))                      # should be 0 NAs
dat <- dat[m1, ]

dat.names <- colnames(dat)                                     # check to make sure desired individuals were obtained 
dat.names
group <- substr(dat.names, 16, 20)
table(group)

y <- DGEList(counts=dat,group=group)                          # Begin formatting for edgeR
dim(y)                                                        # total number of contigs and sample
y$samples                                                     # total number of reads per sample
levels(y$samples$group)                                       # check that all groups are represented

filter <- rowSums(cpm(y)> .5) >= 10                           # keep genes with at least x count per million reads in at least y samples
y.filtered <- y[filter,]
dim(y.filtered)
y.normalized <- calcNormFactors(y.filtered)                   # normalize

plotMDS(y.normalized, labels=group)
plotMDS(y.normalized, labels=group, top = 100)
plotMDS(y.normalized, labels=dat.names, top = 100)            # MDS plot: label options group for group, dat.names for individual sample names

point <- plotMDS(y.normalized, labels=group)

h.group <- which(group == "HFxHM")
w.group <- which(group == "WFxWM")

xs.hatch <- point$x[h.group]
ys.hatch <- point$y[h.group]

xs.wild <- point$x[w.group]
ys.wild <- point$y[w.group]

pdf("Fig1.pdf", onefile=FALSE, width = 7, height=10, paper="letter", title="Figure 4", pointsize=10)
#par(omi=c(4,.1,2,.1), mai=c(1,1,0.5,0.1) ,cex=1.5, las=1)   # here I played with omi and pdf height to raise figure  above center

par(mai=c(5,1.5,1.5,1.5), cex=1.25, cex.axis = 1.05, cex.lab = 1.05, lwd=2)


plot(-10, -10, xlim = c(-0.7, 0.7), ylim = c(-0.6, 0.6), xlab = "Leading Log Fold Change Dimension 1", ylab = "Leading Log Fold Change Dimension 2")
abline(h=0, lwd = 0.5)
abline(v=0, lwd = 0.5)

key=c("Hatch x Hatch", "Wild x Wild" )
legend(0.251,0.634, key, pch=c(16,16), col=c(rgb(230,97,1, maxColorValue=255), rgb(94,60,153, maxColorValue=255)), cex=0.85, bty = "y", box.col="gray65")

points(xs.hatch, ys.hatch, pch=21, bg=rgb(230,97,1, maxColorValue=255), cex=1.75, lwd=0.1)
points(xs.wild, ys.wild, pch=21, bg=rgb(94,60,153, maxColorValue=255), cex=1.75, lwd=0.1)

dev.off()





