# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Performs differetial gene expression usinge edgeR 
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("/home/miles/RNASeq/")
setwd("C:/Dropbox/InPrep/RNASeq/working/DifferentialExpression_MainEffects")
dat.all <- read.csv("genecounts.csv", sep = '\t', header=TRUE)


#======================================================================================================================#


locus <- unlist(strsplit(as.character(dat.all[, 1]),"_"))
loci <- locus[seq(from = 2, to = length(locus), by = 8)]

# Take contig with highest number of hits ===============================================================================#

dat.all2 <- cbind(loci, dat.all)
sumit <- rowSums(dat.all2[, -(1:3)])
dat.all2 <- cbind(dat.all2, sumit)
loc.names <- unique(loci)
counts <- dat.all2[order(dat.all2[, 1], dat.all2[, 93], decreasing = TRUE), ]
counts2 <- counts[match(loc.names, counts[, 1]), ]
counts3 <- counts2[, -c(1,93)]
length(counts3[, 1])

write.table(counts3,file="genecounts.consolidated.mosthits.csv",col.names=TRUE, row.names=FALSE, sep="\t",append=FALSE)


# Take sum of counts across loci =========================================================================================#

dat.all2 <- cbind(loci, dat.all)
loc.names <- unique(loci)
OUT = NULL
for (i in loc.names) {
  counts <- dat.all2[dat.all2[,1] == i, ]
  counts2 <- counts[, -c(1,2)]
  counts3 <- apply(counts2, 2, sum)                   # consider taking the mean
  output <- c(as.character(counts[1,2]), counts3)
  OUT <- rbind(OUT, output)
}
  

write.table(OUT,file="genecounts.consolidated.csv",col.names=TRUE, row.names=FALSE, sep="\t",append=FALSE)