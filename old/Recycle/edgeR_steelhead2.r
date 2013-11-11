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
directory <- "C:/Dropbox/Papers/402_RNASeq/data" 
input <- "GeneCounts.txt"
output <- "DE_genes.txt"


#==================================================================================================================#
setwd(directory)
dat <- read.table(input, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE, row.names=1)

group <- factor(c("hatch","hatch","hatch","hatch","hatch","hatch","hatch","hatch", "hatch","hatch","hatch","hatch",
               "wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild"))               # convert columns to factors


# begin analyses
y <- DGEList(counts = dat, group = group)                                       # create list
y                                                                               # double check list
dim(y)                                                                          # number of "contigs" and samples
y$samples                                                                       # number of contigs ("genes") per sample, aka "tags"

# is the number of unique tags greater than the total number of reads in eachlibrary? if so consider filtering"
# filter           
keep <- rowSums(cpm(y)>1) >= 12  # begin filtering: first number is minimum number of counts per million,
                                 # if 6 million reads then minimum count will be 6 ; second number is across,
                                 # number of samples -usually min group size
y <- y[keep,]
dim(y)

# recompute library sizes
y$samples$lib.size <- colSums(y$counts)
y$samples

#normalize
y <- calcNormFactors(y)
y$samples

#heatmap
#y2 <- predFC(y, prior.count.total=2*ncol(y))  #output can be used for heatmap

#plot MDS of treatments
plotMDS(y)

#s=plotMDS(y)                    #get coordinates for MDS
#asp=cbind(s$x,s$y)
#write.table(asp,file="Males_MDScoord.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)


#estimate dispersion
#prior=getPriorN(y)     #reocmmended in Molec Ecol paper, not sure how to implement

y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)                      #uses empirical Bayes method, if prior null then df=20
#y <- estimateTagwiseDisp(y,prior.df=5)          #I am not sure if this is the correct implementaion, put more weight on tagwise vs. common
y
plotBCV(y)

#test for DE
et <- exactTest(y)
top <- topTags(et)
top
cpm(y)[rownames(top), ]

#total number of genes at 5% FDR
summary(de <- decideTestsDGE(et))    #-1 equals upregulate, 1 equals down regulated:  need to be carefuul though as order (control,vs treatment matters)
top <- topTags(et,n=40) #take top n genes
write.table(top,file="DEgenesboth2.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
