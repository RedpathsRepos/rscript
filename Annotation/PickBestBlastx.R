#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Parses tehe tablular output of hits to Uniprot and reports the highest e-value
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Dropbox/InPrep/RNASeq/working/Annotation/blastx_output")

#======================================================================================================================#

#dat <- read.table(paste(getwd(), "/sprot_output.txt", sep = ''), header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
dat <- read.table(paste(getwd(), "/TrEMBL_output.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)


colnames(dat) <- c("query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit score")

# number of contigs with hits
length(unique(dat[, 1]))
data.frame(table(dat[, 1]))                    # why do some have more than 20?

dat2 <- dat[, c(1, 2, 11)]
dat2 <- dat2[order(dat2[, 1], dat2[, 3]), ]    # just to make sure that ordering is correct

contigs <- unique(dat2[, 1])
m1 <- match(contigs, dat2[, 1])                # remember that match finds the first position and is really fast!

Output <- dat2[m1, ]                            # could easily get top x hits by adusting m1 to have m1 + next x in sequence
write.table(Output, paste(getwd(), "/Top_TrEMBLHits_Contigs.txt", sep = ''), col.names = TRUE, sep="\t", append = FALSE)


# tried all the below first, was overkill and inefficient
#dat3 <- cbind(dat2, 1:length(dat2[, 1]))
#dat4 <- paste(dat3[, 3], "/", dat3[, 4])
#dat4 <- cbind(dat2, dat4)


#by(dat3[, 3:4], dat3[, 1], ReturnMin)   #note by would work well if wanted to take mean of multiple variables (in coloumns) at the same time
#OUT = tapply(dat4[, 4], dat4[, 1], ReturnMin)


#contigs <- unique(dat[, 1])
#wl <- list()
#for (c in contigs){
#wl[[c]] <- dat2[dat2[, 1] == c, ]
#  }





