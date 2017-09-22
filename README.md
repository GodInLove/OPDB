# OPDB
onelink-SRA2OP-DB

## 1.onekey-SRA2OP

<strong>BASIC</strong> V0.0
<hr>
FROM<br/>
`SRA(XXXX)`[1]<br/>
TO<br/>
`fastq-dump (download) (split files)`[6]<br/>
TO[2]<br/>
`.fastq`<br/>
TO[3][4]<br/>
`rockhopper` | `CONDOP` | `RNAseg` <br/>
TO[5]<br/>
`operon_result`<br/>
<hr>
<strong>Extensions:</strong>
[1]considering different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more<br/>
[2]considering the quality of the SRA data<br/>
[3]considering the parameter settings of the software<br/>
[4]considering opting a software artifically or automatically<br/>
[5]considering visual programming technology<br/>
[6]using "split-files" with paired-end data and not using it with single-end data<br/>
