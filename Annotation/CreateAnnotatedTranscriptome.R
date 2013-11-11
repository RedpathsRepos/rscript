#============================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 
# This script:  Merges NR, SPROT, and TrEMBL and blast2go contigs  
# Usage notes:  
#============================================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Dropbox/InPrep/RNASeq/working/Annotation")
library(biomaRt)    # loads bioconductor package

dat.t  <- read.table(paste(getwd(), "/Top_TrEMBLHits_Contigs.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
dat.s  <- read.table(paste(getwd(), "/Top_SprotHits_Contigs.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
dat.n  <- read.table(paste(getwd(), "/Top_NRHits_Contigs.txt", sep = ''), header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
dat.go <- read.table(paste(getwd(), "/b2go_output.txt", sep = ''), header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
contigs<- read.table(paste(getwd(), "/contig.names.txt", sep = ''), header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE) 

#============================================================================================================================#
# Annotate preferentially by Sprot then TrEMBL then NR.  This first step identifies which contigs get annotated by what database.

m1 <- match(dat.t[, 1], dat.s[, 1])                # NA's represent contigs in Trembl that are not in SProt ; isolate for later
dat.t2 <- dat.t[which(is.na(m1)==TRUE), ]

m2 <- match(dat.n[, 1], dat.s[, 1])                # NA's represent contigs in NR that are not in Sprot
dat.n2 <- dat.n[which(is.na(m2)==TRUE), ]
m3 <- match(dat.n2[, 1], dat.t2[, 1])              # NA's represent contigs in NR that are not also in reduced TrEMBL
dat.n2 <- dat.n2[which(is.na(m3)==TRUE), ]

# use bioMart to annotate top SPROT hits=====================================================================================#
vals <- substring(dat.s[, 2], 4, 9)                # Are the uniprot/swissprot acession numbers 
listMarts()                                        # list the avaialble marts (use names from first list, versions in second list)
uniprot <- useMart("unimart")                      # select unimart == uniprot
listDatasets(uniprot)                              # look at the avaialble data sets (only 1 in this case)
uniprot <- useDataset("uniprot", mart = uniprot)   # specifiy the data set you want to use
listFilters(uniprot)                               # shows available filters for the data set
listAttributes(uniprot)                            # attributes define the values that we are interested to retrive
# get gene Id's and basic info
descriptions <- getBM(attributes = listAttributes(uniprot)[c(1:7), 1], filters =listFilters(uniprot)[3, 1] , values = vals , mart = uniprot)  # Do not add more columns of output - risk losing info

# how well did we do?
length(unique(vals))
length(descriptions[, 1])                           # if all ids found a match, these lengths should be equal

# match output back to contig names
head(descriptions)
keep <- cbind(dat.s, vals)
m1 <- match(keep[, 4], descriptions[, 1])
which(is.na(m1) == TRUE)                            # should be 0 
annotated1 <- cbind(keep, descriptions[m1, ])       # annotation with SPROT hits

# use TrEMBL to annotate any extra contigs that were not annotated above=====================================================#
vals <- substring(dat.t2[, 2], 4, 9)
listMarts()                                        # list the avaialble marts (use names from first list, versions in second list)
uniprot <- useMart("unimart")                      # select unimart == uniprot
listDatasets(uniprot)                              # look at the avaialble data sets (only 1 in this case)
uniprot <- useDataset("uniprot", mart = uniprot)   # specifiy the data set you want to use
listFilters(uniprot)                               # shows available filters for the data set
listAttributes(uniprot)                            # attributes define the values that we are interested to retrive
# get gene Id's and basic info
descriptions <- getBM(attributes = listAttributes(uniprot)[c(1:7), 1], filters =listFilters(uniprot)[3, 1] , values = vals , mart = uniprot)  # Do not add more columns of output - risk losing info

# how well did we do?
length(unique(vals))
length(descriptions[, 1])                           # if all ids found a match, these lengths should be equal

# match output back to contig names
head(descriptions)
keep <- cbind(dat.t2, vals)
m1 <- match(keep[, 4], descriptions[, 1])
which(is.na(m1) == TRUE)                            # should be 0 
annotated2 <- cbind(keep, descriptions[m1, ])       # annotation with TrEMBL hits

# use NR to annotate any extra contigs that were not annotated above==========================================================#
# first need to convert gene IDs to Uniprot/Swissport
# couldn't ever get this to work properly
#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#listFilters(mart)
#listAttributes(mart)
#vals <- dat.n2[, 2]
#vals2 <-unlist(strsplit(as.character(vals), '\\|'))
#vals2 <- vals2[seq(from=4, to=length(vals2), 4)]

#results <- getBM(attributes = c("protein_id"), filters = "protein_id", values = vals2, mart = mart)
#results <- getBM(attributes = c("refseq_mrna", "uniprot_swissprot"), filters = "refseq_mrna", values = vals, mart = mart)

# add in blast2go============================================================================================================#
annotated <- rbind(annotated1, annotated2)

# format blast2go - can skip this step if blast2go_formatted.txt already exists; can take several hours
genes <- unique(dat.go[, 1])


OUT = NULL
plot(0,1,xlim=c(0,68361), cex=3)
for (i in 1:length(genes)) {
  gene <- genes[i]
  next.gene <- genes[i+1]
  m1 <- match(gene, dat.go[, 1])
  m2 <- match(next.gene, dat.go[, 1])
  mrows <- m1:(m2-1)
  gene <- dat.go[mrows, ]
  go.ids = paste(unique(gene[,2]), collapse = ' ')
  go.terms = paste(unique(gene[,3]), collapse = ' ')
  output = cbind(as.character(gene[1,1]), go.terms, go.ids)
  OUT = rbind(OUT, output)
  points(i,1,xlim=c(0,68361), cex=3)
}

write.table(OUT, paste(getwd(), "/blast2go_formatted.txt", sep = ''), col.names = TRUE, sep="\t", append = FALSE)

# match annotated to GO output=====================================================================================================#
head(annotated)
head(OUT)
# add 'empty' contigs to go annotation in order to facilitate merging
m1 <- match(annotated[, 1], OUT[, 1])
mis.names <- annotated[which(is.na(m1)==TRUE), 1]
mis.names <- cbind(as.character(mis.names), "NA", "NA")
colnames(OUT) <- c("locus.name", "go.terms", "go.ids")
colnames(mis.names) <- c("locus.name", "go.terms", "go.ids")
OUT2 <- rbind(OUT, mis.names)

m1 <- match(annotated[, 1], OUT2[, 1])
length(which(is.na(m1)==TRUE))     # should be length 0
final1 <- cbind(annotated, OUT2[m1, ])

m1 <- match(OUT2[, 1], annotated[, 1])
m2 <- OUT2[which(is.na(m1)==TRUE),]           # these go-annotated contigs do not match to annotated (match them back to empty contigs and add in)
add.in <- cbind(m2[, 1], "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", m2)
colnames(add.in) <- colnames(final1)
final2 <- rbind(final1, add.in)

write.table(final2, paste(getwd(), "/transcriptome.txt", sep = ''), col.names = TRUE, sep="\t", append = FALSE)

  

# add in empty contigs  - Not used currently
contigs2 <- unlist(strsplit(as.character(contigs[, 1]), ">"))
contigs3 <- contigs2[seq(from = 2, to = length(contigs2), by = 2)]     # all contig names from transcriptome ; could use this to add in Unannotated if wanted


