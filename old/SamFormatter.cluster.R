#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  This script parses a SAM file with a header and creates a new file with (1) no header, (2) the desired columns, (3) rows only for reads that were aligned  
# Usage notes:  Run as is
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("~/Dropbox/Papers/402_RNASeq")
input  <- commandArgs(TRUE)
out.filename <- paste(input,".txt", sep ="")


#======================================================================================================================#

data <- file(input[1], open = "r")
OUT <- NULL

while(length(one.line <- readLines(data, n = 1)) > 0){
  sp.line <- strsplit(one.line, "\t")
  if (length(sp.line[[1]]) < 12) {next} else {     # 1. ignore header lines (lines with less than 12 items)
    out.line <- sp.line[[1]][c(1,2,3,5)]           # 2. extract desired columns
    if (out.line[2] == "4") {next} else {          # 3. ignore reads that were not aligned (designated by "4")
      out.line <- data.frame(t(out.line))
      OUT <<- rbind(OUT, out.line)
    }
  }  
}

close(data)
write.table(OUT, out.filename, col.names = FALSE, sep="\t", append = TRUE) 

