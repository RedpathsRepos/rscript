#==================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script:  Uses edgeR to test for differential gene expression
# Usage notes:  NA
#==================================================================================================================#
# Source files, import packages, set working directory, initialize variables
# source("http://bioconductor.org/biocLite.R")                                  # run these two lines to install package from bioconductor, if not installed
# biocLite("edgeR")
library(edgeR)
directory <- "C:/Dropbox/Papers/402_RNASeq/output" 
input <- "GeneCounts_bowtie.hybrids.txt"
output <- "DE_genes.txt"



#==================================================================================================================#
setwd(directory)
dat <- read.table(input, header=FALSE, sep=" ", na.strings="?", dec=".", strip.white=TRUE, row.names=1)
dat <- dat[-1, ]      #only use this line if there are column headers in the file
dat[, 1:42] <- sapply(dat[, 1:42], as.character) 
dat[, 1:42] <- sapply(dat[, 1:42], as.numeric) #convert to numeric


colnames(dat) <- c("HH1", "HH2", "HH3", "HH4", "HH5", "HH6", "HH7", "HH8", "HH9", "HH10", "HH11", "HW1", "HW2",
                   "HW3", "HW4", "HW5", "HW6", "HW7", "HW8", "HW9", "HW10", "WH1", "WH2", "WH3", "WH4", "WH5",
                   "WH6", "WH7", "WH8", "WH9", "WH10", "WW1", "WW2", "WW3", "WW4", "WW5", "WW6", "WW7", "WW8",
                   "WW9", "WW10", "WW11")   #female designation first (e.g., HW equals Hatch mother and wild father)




group <- c("HH", "HH", "HH", "HH", "HH", "HH", "HH", "HH", "HH", "HH", "HH", "HW", "HW", "HW", "HW", "HW"
                  , "HW", "HW", "HW", "HW", "HW", "WH", "WH", "WH", "WH", "WH", "WH", "WH", "WH", "WH", "WH"
                  , "WW", "WW", "WW", "WW", "WW", "WW", "WW", "WW", "WW", "WW", "WW")               # neccesssary grouping information



# Begin formatting for edgeR
y <- DGEList(counts=dat,group=group)
dim(y)  # total number of contigs and sample
y$samples  # total number of reads per sample
levels(y$samples$group)   # check that all groups are represented


#filtering
filter <- rowSums(cpm(y)>1) >= 4 #kepp genes with at least 1 count per million reads in at least 4 samples
y.filtered <- y[filter,]
dim(y.filtered)

#normalize
y.normalized <- calcNormFactors(y.filtered)

#data exploration
barplot(y.normalized$samples$lib.size*1e-6, ylab = "library size (millions)")  #creates barplot of filtered read #
plotMDS(y.normalized)
s=plotMDS(y.normalized)                    #get coordinates for MDS
asp=cbind(s$x,s$y)
write.table(asp,file="hybrids_all_MDScoord.txt",col.names=TRUE, row.names=TRUE, sep="\t",append=FALSE)



counts <- cpm(y.normalized, normalized.lib.sizes=TRUE)
pc <- prcomp(counts)
plot(pc$rotation[1:11, 1], pc$rotation[1:11, 2], pch=21, bg="orange")
points(pc$rotation[12:21, 1], pc$rotation[12:21, 2], pch=21, bg="green")
points(pc$rotation[22:31, 1], pc$rotation[22:31, 2], pch=21, bg="yellow")
points(pc$rotation[32:42, 1], pc$rotation[32:42, 2], pch=21, bg="blue")

heatmap(t(counts[1:100,]),Rowv=NA,Colv=NA,col = heat.colors(160), scale="column", margins=c(5,10),xlab="Gene Ids",ylab="Fish Type")
# using the glm approach, requires a design matrix
# read in and format the design matrix
design <- read.table("HybridDesignMatrix.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
matrix <- factor(design[, 2])
treatment <- factor(design[, 3])
design <- model.matrix(~0+matrix+treatment)
rownames(design) <- colnames(y.normalized)

# estimate the dispersion
y.processed <- estimateGLMCommonDisp(y.normalized, design, verbose = TRUE)  # estimate the overall dispersion
y.processed <- estimateGLMTrendedDisp(y.processed, design)    # estimate gene-wise dispersion rates
y.processed <- estimateGLMTagwiseDisp(y.processed, design)
plotBCV(y.processed)

# fit the genewise glms
fit <- glmFit(y.processed, design)

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
#lrt <- glmLRT(fit)
topTags(lrt)
FDR <- p.adjust(lrt$table$PValue, method = "BH")
sum(FDR < 0.05)

top <- rownames(topTags(lrt, 150))
asp <- cpm(y.processed)[top,]
summary(dt <- decideTestsDGE(lrt))

isDE <- as.logical(dt)
de.names <- rownames(y.processed)[isDE]
de.names  # names of the DE genes

plotSmear(lrt, de.tags = de.names)
abline(h=c(-1,1), col = "blue")   #lines equal 2 fold up or down

