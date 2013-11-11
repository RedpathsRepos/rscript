#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  This script parses a SAM file with a header and creates a new file with (1) no header, (2) the desired columns, (3) rows only for reads that were aligned  
# Usage notes:  Run as is, the fill, quote, and colClasses are all neccessary in the read.table.call
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("~/Dropbox/Papers/402_RNASeq")
library("R.utils")
header.stop <- 130602  #first line where data starts
chunk.size <- 100000   #chunk size in number of lines
col.classes <- c("character", "character", "character", "NULL", "character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL")
#======================================================================================================================#

nlines <- countLines(paste(getwd(), "/Data/test.sam", sep = ''))
start <- seq(from = header.stop, to = nlines, by = chunk.size)
chunks <- c(rep(chunk.size,length(start)-1), nlines-start[length(start)])

OUT = NULL
OUT2 = NULL
for (p in 1:length(start)){
  data <- read.table(paste(getwd(), "/Data/test.sam", sep = ''), header=FALSE, sep="\t", na.strings="NA", dec=".", fill = TRUE, quote = "", strip.white=TRUE, skip = start[p], nrows = chunks[p], colClasses=col.classes)  
  output <- data[data[, 2] != 4, ]
  OUT = rbind(OUT,output)
}

write.table(OUT, paste(getwd(), "/output/output.txt", sep = ''), col.names = FALSE, sep="\t", append = TRUE) 


