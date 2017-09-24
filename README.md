### Usage

[download the whole package] and decompress it

test: 
	python SRA2OP_0_0.py SRR5486953

**BASIC** V0.0
___
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

[download the whole package]:https://github.com/GodInLove/OPDB.git
