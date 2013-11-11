setwd("C:/POPS/GeneExpression")
dat<- read.table("blastx_conservative.txt", header=FALSE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)

#grep(dat[,1]=="<Hit_def>")
#col1=grep("<Hit_def>",dat[,1])
#col1=dat[col1,]

OUT=NULL
col1=grep("<BlastOutput_query-def>",dat[,1])
for (i in col1) {
 data=dat[i:(i+35),]
 col2=data[grep("<Hit_def>",data)]
 col3=data[grep("<Hsp_evalue>",data)]
 out=cbind(as.character(data[1]),as.character(col2),as.character(col3))
 if(length(out)==1) {out=cbind(as.character(data[1]),"null","null")}
 OUT=rbind(OUT,out)
 }
 
write.table(OUT,file="blast_results_females.txt",col.names=TRUE,sep="\t",append=FALSE)