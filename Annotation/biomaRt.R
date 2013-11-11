#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Takes a blastx output file (table output format) and usese biomart for annotation purposes  
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Dropbox/InPrep/RNASeq/working/Annotation")


library(biomaRt)    # loads bioconductor package
#source(paste(getwd(), "/source/.R", sep = ''))


#======================================================================================================================#

#dat <- read.table(paste(getwd(), "/Top_SprotHits_Contigs.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
dat <- read.table(paste(getwd(), "/Top_TrEMBLHits_Contigs.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
values <- substring(dat[, 2], 4, 9)

listMarts()                     # list the avaialble marts (use names from first list, versions in second list)
uniprot <- useMart("unimart")   # select unimart == uniprot

listDatasets(uniprot)           # look at the avaialble data sets (only 1 in this case)
uniprot <- useDataset("uniprot", mart = uniprot)  # specifiy the data set you want to use

listFilters(uniprot)            # shows available filters for the data set
listAttributes(uniprot)         # attributes define the values that we are interested to retrive

vals = values
# get gene Id's and basic info
descriptions2 <- getBM(attributes = listAttributes(uniprot)[c(1:7), 1], filters =listFilters(uniprot)[3, 1] , values = vals , mart = uniprot)  # note that there are ~33,000 hits because multiple contigs had the exaxt same protein acession number (can switch this option off if needed: query ?getBM)

# match output back to contig names
head(descriptions2)
keep <- cbind(dat, values)

m1 <- match(keep[, 4], descriptions2[, 1])
which(is.na(m1) == TRUE)         # should be 0 for descriptions2
annotated1 <- cbind(keep, descriptions2[m1, ])


# find and add in GO terms, if avaialble
descriptions1 <- getBM(attributes = listAttributes(uniprot)[c(1,8,9,10,12), 1], filters =listFilters(uniprot)[3, 1] , values = vals , mart = uniprot)  # more rows, but fewer hits...typically

# the below may take a long time, but it only took about 10 mins for 6000 GO values

annotated2 <- cbind(annotated1, 0)
go.values <- unique(descriptions1[, 1])
length(go.values)
for (i in go.values){
  go1 <- descriptions1[descriptions1[, 1]==i, 2]
  go2 <- paste(go1, collapse = ',')     
  find.locations <- which(annotated2[, 4] == i)
  annotated2[find.locations, 12] <- go2   # where 12 == ncol(annotated2)
}


write.table(annotated2, paste(getwd(), "/TrEMBL_annotations.txt", sep = ''), col.names = TRUE, sep="\t", append = FALSE)
