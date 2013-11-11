setwd("C:/POPS/GeneExpression")
dat<- read.table("blastx_moderate.txt", header=FALSE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)

#grep(dat[,1]=="<Hit_def>")
#col1=grep("<Hit_def>",dat[,1])
#col1=dat[col1,]


col1=grep("<BlastOutput_query-def>",dat[,1])
col2=grep("<Parameters_expect>",dat[,1])
#col3=grep("<Parameters_gap",dat[,1])
col4=grep("<Hit_len>",dat[,1])
col5=grep("<Hit_num>",dat[,1])
col6=grep("<Hit_id>",dat[,1])
col7=grep("<Hit_def>",dat[,1])

ids=sort(c(col1,col2,col4,col5,col6,col7))

dats=dat[ids,]

OUT=NULL
for (i in 1:length(col1)) {
 starts=col1[i]
 stops=col1[i+1]
 id1=ids[which(ids>=starts)]
 id2=id1[which(id1<stops)]

 if(length(id2)==6) {out=t(as.character(dat[id2,1]))
    out=cbind(as.character(substr(out[,1],24,31)),as.character(substr(out[,2],20,24)),as.character(substr(out[,3],10,10)),as.character(out[,4]),as.character(out[,5]),as.character(substr(out[,6],9,100)))
 }
 if((length(id2)<6)==TRUE) {
    out=t(dat[id2,1])
    out=cbind(as.character(substr(out[,1],24,31)),as.character(substr(out[,2],20,24)),"","","","")}

 OUT=rbind(OUT,out)
 }      #should get error in NextMethod...

 
colnames(OUT)<-c("Contig","Evalue","Nhits","Hit Id","Hit def","Hit Length")

dat<- read.table("DEgenesFemales.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)

m1=match(OUT[,1],as.character(dat[,1]))
OUT2=cbind(dat[m1,],OUT)

write.table(OUT2,file="Blastx_results_females_moderate.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)