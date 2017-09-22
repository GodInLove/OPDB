# OPDB
onelink-SRA2OP-DB

## 1.onekey-SRA2OP

<strong>BASIC</strong> V0.0
<hr>
FROM  
`SRA(XXXX)`[1]
TO
`fastq-dump (download) (split files)`[6]
TO[2]
`.fastq`
TO[3][4]
`rockhopper` | `CONDOP` | `RNAseg`
TO[5]
`operon_result`
<hr>
<strong>Extensions:</strong>
[1]considering different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
[2]considering the quality of the SRA data
[3]considering the parameter settings of the software 
[4]considering opting a software artifically or automatically
[5]considering visual programming technology
[6]using "split-files" with paired-end data and not using it with single-end data
