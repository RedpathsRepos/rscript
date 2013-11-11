#==================================================================================================================#
# Script created by Mark Christie, contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script: is used to create an index for Rsubread ; 
# Usage notes:  installed on CGRB cluster only, run as SGE_Batch -c 'myR --vanilla < "script.r"> ./log.txt' -o StdOut
# Usage notes:  Also make sure to run on locally installed R; bash myR
#==================================================================================================================#
# Source files, import packages, set working directory, initialize variables
transcriptome <- "lemma/christie/utlities/steelhead_transcriptome_v1.86402.fa"
output.dir <- "/lemma/christie/subread/index" 

#==================================================================================================================#
setwd(output.dir)
library("Rsubread", lib="~/local_programs/R-3.0.1/packages/")

buildindex(basename = "./salmon_index", reference = transcriptome)

