#======================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Matches annotated transcriptome back to edgeR results  
# Usage notes:  run line by line
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Dropbox/InPrep/RNASeq/working")


#library(biomaRt)    # loads bioconductor package
#source(paste(getwd(), "/source/.R", sep = ''))


#======================================================================================================================#

annotation <- read.table(paste(getwd(), "/uniprot_annotations.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
de.genes <- read.table(paste(getwd(), "/topgenes_new.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
de.gene.names <- rownames(de.genes)   # this step is needed because de.genes was saved with rownames=TRUE

m1 <- match(de.gene.names, annotation[, 1])
id <- cbind(de.genes, annotation[m1, ])

de.gene.names2 <- substring(de.gene.names, 1,12)
annotations2 <- substring(annotation[, 1], 1,12)

m1 <- match(de.gene.names2, annotations2)
id <- cbind(de.genes, annotation[m1, ])

sort(id$protein_name)

go.id <- id[, 17]
go.id1 <- go.id[-(which(go.id==0))]   # remove 0's, those with a gene name but no GO annotation
go.id2 <- go.id1[-(which(is.na(go.id1)==TRUE))]
ngo <- length(go.id2)
OUT <- NULL
for (i in 1:ngo){
  set <- go.id2[i]
  set2 <- strsplit(as.character(set), ',')
  set3 <- as.data.frame(set2)
  set3 <- cbind(i, set3)
  colnames(set3) <- c("id", "GO")
  OUT <- rbind(OUT, set3)
}

OUT <- unique(OUT[, 2])
write.table(OUT, paste(getwd(), "/GO_Ids.txt", sep = ''), col.names = TRUE, sep="\t", append = FALSE)





