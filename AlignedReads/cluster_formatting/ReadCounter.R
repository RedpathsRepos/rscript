#==================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script:  Takes a sam file and counts all the reads 
# Usage notes:  column 2 is the SAM Flag, reads with 0 are mapped to the forward strand, reads with 16 are mapped
# to the reverse strand, and reads with 4 are unmapped, other values can be meaningful so search for those as well.
#==================================================================================================================#
# Source files, import packages, set working directory, initialize variables
# consider adding chunk_Reader_regularfile  to work with large files
directory <- "C:/Dropbox/Papers/402_RNASeq/R_scripts" 
input1 <- "contigs.txt"
input2 <- "hybrid07.fastq.sam.txt"
output <- "counts.txt"


#==================================================================================================================#
setwd(directory)
contigs <- read.table(input1, header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
counts <- read.table(input2, header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)  #could result in memory issues if very lager

counted <- data.frame(table(counts[,2]))

match.contigs <- match(contigs[,1],counted[,1])
zero.counts <- which(is.na(match.contigs)=="TRUE")
zero.counts2 <- cbind(as.character(contigs[zero.counts,]), 0)

colnames(zero.counts2) <- c("Contig", "Count")
colnames(counted) <- c("Contig", "Count")
counted <- rbind(counted, zero.counts2)

match.to.contigs <- match (contigs[,1], counted[,1])  #preserve the exact same row order as in the contigs file
counted <- counted[match.to.contigs, ]

write.table(counted, output, col.names = TRUE, row.names = FALSE, sep="\t", append = FALSE)  



