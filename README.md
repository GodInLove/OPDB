___
### Usage
[download the whole package] and decompress it

**test:** 
	
	python SRA2OP_0_0.py SRR5486953

**dependencies:**
1. Java version 1.6 or later
2. Python 2.7
___
**BASIC** V0.0
FROM

`SRA(XXXX)`*1*

TO

`fastq-dump(download)(split files)`*2*

TO

`.fastq`

TO *3* *4*

`rockhopper` | `CONDOP` | `RNAseg`

TO

`operon_result`
___
**Extensions:**
1. considering different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
2. considering the quality of the SRA data
3. considering the parameter settings of the software
4. considering opting a software artifically or automatically
5. considering visual programming technology
6. using "split-files" with paired-end data and not using it with single-end data
___
**Still cannot resolve it**
1. the annotion file(`.gff`|`.rnt`|`.ptt`|`.fna`|`.opr`)
2. the Rscript(How it can work like the bash-shell?For example, `rm *.fna`)
3. It will not go on to finish its job if there are some bugs,although I have fix some bugs. 
4. It still cannot automaticlly choose the software(`rockhopper`|`CONDOP`|`RNAseg`)
5. It still cannot filter the bad SRR data

[download the whole package]:https://github.com/GodInLove/OPDB.git
