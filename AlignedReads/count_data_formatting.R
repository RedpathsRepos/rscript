# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Takes processed Sam files, merges them onto a set of desired contigs, orders individuals by column, and prints desired output file  
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
# setwd("/home/miles/RNASeq/")
setwd("C:/Dropbox/Papers/402_RNASeq/working")
#library()
source("C:/Dropbox/Rscript/FileImporter.R")
data.directory <- paste(getwd(), "/alignedreads/", sep = '')

#======================================================================================================================#

all.data <- FileImporter(directory = data.directory)           # import all files into a single list
names(all.data)                                                # see the names of each list element

names <- sort(all.data[[1]][, 1])                              # sort the locus names 
new.data <- matrix(nrow = length(names),ncol = 91)             # create an empty matrix
new.data[, 1] <- as.character(names)                           # add locus names to first column
 
for (i in 2:length(all.data)) {                                # start at 2 becuase the very first item in the list is simply the contig names
  m1 <- match(all.data[[1]][, 1], all.data[[i]][, 1])          # for list item i, match to locus name
  zero.reads <- which(is.na(m1) == TRUE)                       # identify which loci had no match and thus had 0 reads
  loc.names <- all.data[[1]][zero.reads, 1]                    # get the names of these loci
  new.loci <- data.frame(loc.names, 0)                         # create new data with locus name and 0, for 0 counts
  colnames(new.loci) <- c("V1", "V2")     
  new.loci2 <- rbind(all.data[[i]], new.loci)                  # attach loci with 0 reads to the other loci
  m2 <- match(new.data[, 1], new.loci2[, 1])                   # match up row (locus) names
  print(length(which(is.na(m2)==TRUE)))                        # serves as a check, should all be 0 as all loci should be included
  new.data[, i] <- new.loci2[m2, 2]                            # add correctly indexed counts to new.data
}

#begin ordering of samples by treatments, matrices, and fry sex



dat <- read.csv("database.names.csv")
dat <- dat[order(dat[, 3], dat[, 5], dat[, 2]), ]                # order by treatment, sex, and then matrix  
m1 <- match(dat[, 1], substr(names(all.data), 1, 14))            # match orderd names to position in new.data, note substr is used to remove .sam.txt

names(all.data)                                                  # give new.data column names, these first two lines the order should be identical, serve as check
c("contig.names", as.character(sort(dat[, 6])))
colnames(new.data) <- c("contig.names", as.character(sort(dat[, 6])))
new.data <- new.data[, c(1,m1)]                                  # give approriate column order to new.data
cbind(colnames(new.data)[-1], as.character(dat[, 6]))            # check that new order matches with wanted order


write.table(new.data, "genecounts.csv", row.names = FALSE, col.names = TRUE, sep="\t", append = FALSE)   


