___
### Usage
[download the whole package] and decompress it

**test:** 
	
	python main.py -i SRR5486953 -o /home/lyd/Desktop -m 1 -k eco
	python main.py -h

**dependencies:**
1. Java version 1.6 or later
2. Python 2.7
___
**BASIC** V0.0
FROM

`SRA(XXXX)`

TO

`fastq-dump(download)(split files)`

TO

`.fastq`

TO

`rockhopper` | `CONDOP` | `RNAseg`

TO

`operon_result`
___
**Extensions:**
1. get the information about different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
2. get the information about layout:PAIRED-END or SINGLE-END, then automaticlly set the pramater. 
3. offer the options to choose the software and processor number.
___
**Still cannot resolve it**
1. the annotion file(`.gff`|`.rnt`|`.ptt`|`.fna`)
2. It will not go on to finish its job if there are some bugs,although I have fix some bugs. 
3. I have not yet tested the software called RNAseg because It took too much memory size.
4. It still cannot filter the bad SRR data
5. It did not filter the fastq file which had bad quality.
6. It is time to consider visual programming technology

[download the whole package]:https://github.com/GodInLove/OPDB.git
