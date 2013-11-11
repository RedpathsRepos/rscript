#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  This script parses a SAM file with a header and creates a new file with (1) no header, (2) the desired columns, (3) rows only for reads that were aligned  
# Usage notes:  Run as is
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("~/Dropbox/Papers/402_RNASeq")

#======================================================================================================================#

data <- file(paste(getwd(), "/Data/test.sam", sep = ''), open = "r")
OUT <- NULL
n.line <- 0
while(length(one.line <- readLines(data, n = 1)) > 0){
  sp.line <- strsplit(one.line, "\t")
  n.line <- n.line+1
  print(n.line)
  if (length(sp.line[[1]]) < 12) {next} else {     # 1. ignore header lines (lines with less than 12 items)
    out.line <- sp.line[[1]][c(1,2,3,5)]           # 2. extract desired columns
    if (out.line[2] == "4") {next} else {          # 3. ignore reads that were not aligned (designated by "4")
      out.line <- data.frame(t(out.line))
      OUT <<- rbind(OUT, out.line)
    }
  }  
}

close(data)
write.table(OUT, paste(getwd(), "/output/output.txt", sep = ''), col.names = FALSE, sep="\t", append = TRUE) 


