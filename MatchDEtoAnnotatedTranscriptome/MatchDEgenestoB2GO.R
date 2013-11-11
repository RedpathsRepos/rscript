# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Matches DE genes to annotated transcriptome or B2GO
# Usage notes:  
#======================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
#setwd("/home/miles/RNASeq/")
setwd("C:/Dropbox/InPrep/RNASeq/working/MatchDEtoAnnotatedTranscriptome")

dat1 <- read.table("TopGenes_AllContigs.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
b2go <- read.table("b2go_output2.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)

#======================================================================================================================#
length(b2go[, 1])
length(unique(b2go[, 1]))

dat1 <- dat1[dat1[, 4] <= 0.05, ]   # filter by p-value?
dat1 <- dat1[1:100, ]


genes <- rownames(dat1)
genes2 <- unlist(strsplit(genes, "T"))
genes2 <- genes2[seq(1, length(genes2), by =2)]

bgenes2 <- unlist(strsplit(as.character(b2go[,1]), "T"))
bgenes2 <- bgenes2[seq(1, length(bgenes2), by =2)]
b2go <- cbind(bgenes2, b2go)


m1 <- match(genes, b2go[, 2])   # match to exact trancript ; # match ony finds first position, need to use for loop
OUT = NULL
for (i in 1:length(genes)){
  gene <- genes[i]
  output <- b2go[b2go[, 2]==gene, ]
  if(is.na(output[1,1])==FALSE) {output <- cbind(output, dat1[i, 5])}
  OUT <- rbind(OUT, output)
}

m1 <- match(genes2, b2go[, 1])  # match to same locus ; # match ony finds first position, need to use for loop
OUT = NULL
for (i in genes2){
  gene <- i
  output <- b2go[b2go[,1]==gene, ]
  OUT <- rbind(OUT, output)
}



write.table(OUT, paste(getwd(),"/GOterms.top100.txt", sep = ''), col.names = FALSE, sep="\t", append = FALSE)

