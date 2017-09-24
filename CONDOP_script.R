Args <- commandArgs(TRUE)
argv <- Args[1]

library(seqinr)
library(mclust)
library(CONDOP)

ct <- read.csv(paste(argv,"_output/",argv,"_count",sep=""),header = TRUE)
file_genome_annot <- paste(argv,"_ref/.gff",sep="")
file_operon_annot <- paste(argv,"_ref/.opr",sep="")
file_fna_annot <- paste(argv,"_ref/.fna",sep="")

data.in <- pre.proc(file_genome_annot,file_operon_annot,file_fna_annot,list.cov.dat = list(ct = ct))

res.comapa <- run.CONDOP(data.in, bkgExprCDS=0.2, bkgExprIGR=0.2,maxLenIGR=150, win.start.trp=c(100,10),win.end.trp=c(10,100), norm.type = "n1",cl.run = 30, nfolds = 5, cons = 2,find.ext=TRUE,return.all = FALSE,save.TAB.file = "CONDOP",verbose = TRUE)