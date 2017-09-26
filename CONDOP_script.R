Args <- commandArgs(TRUE)
srr_n <- Args[1]
ref <- Args[2]
output <- Args[3]

library(seqinr)
library(mclust)
library(CONDOP)
library(stringr)

filename <- dir(ref)
fna <- grep(".fna",filename,TRUE)
gff <- grep(".gff",filename,TRUE)
opr <- grep(".opr",filename,TRUE)

ct <- read.csv(paste(output,"/",srr_n,"_count",sep=""),header = TRUE)
file_genome_annot <- paste(ref,"/",gff,sep="")
file_operon_annot <- paste(ref,"/",opr,sep="")
file_fna_annot <- paste(ref,"/",fna,sep="")

data.in <- pre.proc(file_genome_annot,file_operon_annot,file_fna_annot,list.cov.dat = list(ct = ct))

res.comapa <- run.CONDOP(data.in, bkgExprCDS=0.2, bkgExprIGR=0.2,maxLenIGR=150, win.start.trp=c(100,10),win.end.trp=c(10,100), norm.type = "n1",cl.run = 30, nfolds = 5, cons = 2,find.ext=TRUE,return.all = FALSE,save.TAB.file = "CONDOP",verbose = TRUE)