#requires a table of gene counts as input (i.e., use bowtie then counter.pl) - Note can only work for 2 treatments currently.

library(NBPSeq)
#data(arab)  #using built in data as example

setwd("C:/POPS")
arab<- read.table("GeneCounts.txt", header=F, sep=" ", na.strings="?", dec=".", strip.white=TRUE,row.names=1)             #Note use of row.names here!! - sets [1,1] to first data value  - try as.matrix if need
arab=as.matrix(arab)


#workflow 1
#grp.ids=c(1,1,1,2,2,2)
#prepped=prepare.nbp(arab,grp.ids,print.level =5)
 #prep2 = estimate.disp(prepped, method ="NBP", print.level =5)
#out = exact.nb.test(


#workflow 2 - default
grp.ids=c(1,1,2,1,1,2,2,1,2,2,1,1,2,1,2,2,2,1,2,1,1,2,1,2)  # max of 2 groups; 1==wildxwild, 2==hatchxhatch
grp1=1
grp2=2

out=nbp.test(arab, grp.ids, grp1, grp2, print.level=5)
attributes(out)

#count the number of differentially expressed genes
alpha = 0.05
sig.out = out$q.values < alpha
table(sig.out)

significant=which(out$q.values < alpha)
qvalue=out$q.values[significant]
expression1=out$expression.levels[significant,]
fold=out$log.fc[significant]

Output=arab[significant,]
Output2=cbind(Output,data.frame(fold),qvalue)
write.table(Output2,file="Qvalues.txt",col.names=TRUE,sep="\t",append=FALSE)  #move column headers one column to the right (inapproporiate format due to as.matrix)
