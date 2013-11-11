#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Parses the multiple tabular output produced from blastx of the steelhead trancriptome _from NR database after running BlastParse.pl on each xml file 
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Dropbox/InPrep/RNASeq/working/parse_blastxoutput")

#library()
#source(paste(getwd(), "/source/.R", sep = ''))

#======================================================================================================================#

#write.table(output, paste(getwd(), "/output/.txt", sep = ''), col.names = FALSE, sep="\t", append = TRUE)   

directory <- "C:/Dropbox/InPrep/RNASeq/working/parse_blastxoutput/"
files <- list.files(path = directory)
print(files)

mycols <- c("factor", "NULL", "NULL", "character", "NULL", "NULL", "NULL", "NULL", "character", "NULL", "NULL", "NULL", "NULL", "NULL")

OUT <- NULL
for (i in 1:length(files)) {
  name <- files[i]  
  data  <- read.table(paste(directory, files[i], sep = ''), header=TRUE, sep="\t", colClasses=mycols ,quote = "",fill = TRUE, na.strings="?", dec=".", strip.white=TRUE) # note use of fill = TRUE ; turned out to be very important here
  OUT <- rbind(OUT, data)
}

dat <- OUT

# number of contigs with hits
length(unique(dat[, 1]))
data.frame(table(dat[, 1]))                    # why do some have more than 20?

dat2 <- dat[order(dat[, 1], dat[, 3]), ]       # order by locus name and then evalue
#dat2 <- dat[, c(1, 3, 4, 9)]

contigs <- unique(dat2[, 1])
m1 <- match(contigs, dat2[, 1])                # remember that match finds the first position and is really fast!

Output <- dat2[m1, ]                            # could easily get top x hits by adusting m1 to have m1 + next x in sequence
write.table(Output, paste(getwd(), "/Top_NRHits_Contigs.txt", sep = ''), col.names = TRUE, sep="\t", append = FALSE)



