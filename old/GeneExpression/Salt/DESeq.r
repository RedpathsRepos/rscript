#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")                                 #run these 2 lines only to install

library("DESeq")

setwd("C:/POPS/Salt")
dat <- read.table("GeneCounts.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE,row.names=1)      #use header with treatment and id (wild1,wild2,etc...)
conds <- factor(c("control","control","control","salt","salt","salt"))

cds <- newCountDataSet(dat, conds)
cds <- estimateSizeFactors(cds)           #can be used to standardize samples, make them relative to one another, divide raw counts by these factors
sizeFactors(cds)
head(counts(cds,normalized=TRUE))          #"normalized" dataset

cds <- estimateDispersions(cds)          #estimate dispersion
str(fitInfo(cds))

#visualize dispersion:
plotDispEsts <- function( cds ) {
  plot(
  rowMeans( counts( cds, normalized=TRUE ) ),
  fitInfo(cds)$perGeneDispEsts,
  pch = '.', log="xy" )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
  }

plotDispEsts(cds)    # plot dispersion

#test for differential gene expression:
res <- nbinomTest( cds, "control", "salt" )        #note that 1 and 2 are used here becuase those are how I coded the treatements:  change to wild and hatch...
head(res)
#plot
plotDE <- function( res )
 plot(
 res$baseMean,
 res$log2FoldChange,
 log="x", pch=20, cex=.3,
 col = ifelse( res$padj < .1, "red", "black" ) )        #the .1 sets the FDR to 0.1 for plotting red points

plotDE( res )

#create histogram
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")      #may expect a peak on the low end if there are lots of significant genes

#filter for most significant genes
resSig <- res[ res$padj < 0.2, ]
head( resSig[ order(resSig$pval), ], 100 )   # all genes (down or up regulated
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )  #most strongly downregulated


write.table(resSig,file="Pvalues.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)


#Now try with filtering

rs = rowSums(counts(cds))
use <- (rs > quantile(rs, 0.4))
table(use)

cds <- newCountDataSet(dat, conds)
cds <- cds[ use, ]     #go back to step x