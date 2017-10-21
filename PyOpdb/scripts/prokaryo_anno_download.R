###Download gff/ptt/rnt annotation files of a prokaryotic genome ####

## Escherichia coli str. K-12 substr. MG1655 (E. coli)
## KEGG ID is 'eco' (http://www.genome.jp/kegg-bin/show_organism?org=eco)
## NCBI assembly ID is 'GCF_000005845.2' (https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2)

Args <- commandArgs(TRUE)

kegg_id <- Args[1]

## temporary folder for saving files
saveFolder <- Args[2]

library('ProGenome') ## version >= 0.06
library('magrittr')


## check folder
if (!dir.exists(saveFolder)) {
  dir.create(saveFolder)
} else {}

## download gff and feature_table files
# download.SpeAnno(kegg_id, 'gff', saveFolder)
# download.SpeAnno(kegg_id, 'feature_table', saveFolder)
# download.SpeAnno(kegg_id, '[^from]_genomic.fna', saveFolder)
files <- dir(saveFolder, full.names = TRUE)


## extract the fna file as 'eco.fna'
# files %>%
#   grepl('[^from]_genomic.fna', .) %>%
#   `[`(files, .) %>%
#   paste0('zcat ', ., ' > ', file.path(saveFolder, paste(kegg_id,".fna",sep=""))) %>%
#   system
## extract the gff file as 'eco.gff'
# files %>%
#   grepl('gff', .) %>%
#   `[`(files, .) %>%
#   paste0('zcat ', ., ' > ', file.path(saveFolder, paste(kegg_id,"_orgin.gff",sep=""))) %>%
#   system

## extract ptt files as 'eco.ptt'
## check ptt at http://cs.wellesley.edu/~btjaden/genomes/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.ptt
## extract rnt files as 'eco.rnt'
## check rnt at http://cs.wellesley.edu/~btjaden/genomes/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.rnt
ft <- files %>%
  grepl('feature', .) %>%
  `[`(files, .) %>%
  read.gff

ft %>%
  ExtractPtt %>%
  write.ptt(file.path(saveFolder, paste(kegg_id,".ptt",sep="")))

ft %>%
  ExtractRnt %>%
  write.rnt(file.path(saveFolder, paste(kegg_id,".rnt",sep="")))

# paste0("rm ",saveFolder,"/*.txt")
# paste0("rm ",saveFolder,"/*.gz")