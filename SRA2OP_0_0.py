#!bin/python2.7
import os
import getopt
import sys
import re
import urllib2

# Basic V0.0

def usage():
	print "Useage:python SRA2OP_0_0.py SRR5486953"

#1.input SRRXXXX and check it
def checkargv(argv):
	pat = re.compile(r'[SRA]+[0-9]+')
	searchobj = re.search(pat,argv)
	if not searchobj:
		print "wrong SRAnumber!"
		usage()
		sys.exit(2)

def mkdir(argv):
	os.system("mkdir "+argv+"_input")
	os.system("mkdir "+argv+"_ref")
	os.system("mkdir "+argv+"_output")

def not_qualified():
	print "the SRR data is not qualified!"
	sys.exit(2)

def check_SRR(information):
	#check whether the SRR data is up to standard or not
	if information['Organism'] == "Escherichia coli":
		not_qualified()
	else:
		pass

#2.figure out something and prepare reference data 
'''
figuring out different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
figuring out the parameter settings of the software
figuring out whether the data is paired-end or not
figuring out .opr file
figuring out .gff .rnt .fna .gtf
'''
#WRITE a python/R/Perl script to connect NCBI API and get those information
def getinformation(argv,information):
	retry = -10
	while retry < 0:
		# print retry
		try:
			url = 'https://www.ncbi.nlm.nih.gov/sra/?term='+argv+'%5BAll+Fields%5D+AND+"biomol+rna"%5BProperties%5D'
			print "open link:\n"+url+"\nplease wait...\n"
			page = urllib2.urlopen(url,timeout=10)
			html = page.read()
			retry = 1
			pat_notfind = re.compile(r'No items found')
			searchobj = re.search(pat_notfind,html)
			if searchobj:
				print "SRA Experiment "+argv+" is not public or not RNA. please change your SRAnumber"
				sys.exit(2)
			#find Organism
			pat_org = re.compile(r'Organism\: \<span\>\<a href\=.*\"\>([a-zA-Z0-9\.\s\-]+).*expand showed sra-full-data')
			org = re.findall(pat_org,html)
			if len(org)==1:
				information['Organism'] = org[0]
			else:
				print "Cannot find the Organism.WHAT happen?the result is",org
				sys.exit(2)
			#find Instrument
			pat_ins = re.compile(r'Instrument\: \<span\>([A-Za-z0-9\s]+)\<\/span\>\<\/div\>\<div\>Strategy\:')
			ins = re.findall(pat_ins,html)
			if len(ins)==1:
				information['Instrument'] = ins[0]
			else:
				print "Cannot find the Instrument.WHAT happen?the result is",ins
				sys.exit(2)
			#find Layout
			pat_lay = re.compile(r'Layout: <span>([PAIREDSINGLEpairedsinglenN]+).*sra-full-data')
			lay = re.findall(pat_lay,html)
			if len(lay)==1:
				information['Layout'] = lay[0]
			else:
				print "Cannot find the Layout.WHAT happen?the result is",org
				sys.exit(2)
			check_SRR(information)
			f = open(argv+"_information.txt",'w')
			# print information
			print "\nNow you can open your workspace, and you will find a file named "+argv+"_information.txt. It stores some BASIC information about your SRR data."
			for item in information:
				f.write(item+": "+information[item]+"\n")
			f.close()
		except:
			retry+=1
			if retry == 0:
				print "wrong!!Cannot connect NCBI"
				sys.exit(2)

#WRITE a python/R/Perl script to connect NCBI API and download the ref.fna ref.gff ref.rnt ref.gtf
def download_annotion(argv,organism):
	print "downloading annotion:.gtf .rnt .gff......."
	# assume that finished it
	# download the annotion in the ref_dir
	os.system("cp NC_000913_ref/* "+argv+"_ref/")
	print organism+"annotion files were downloaded in the "+argv+"_ref"

#WRITE a python/R/Perl script to figure out whether it is paired-end or not
def paired_or_single(layout):
	pat_paired = re.compile(r'[padPAD]')
	searchobj = re.search(pat_paired,layout)
	if searchobj:
		return 1
	else:
		return 0

#WRITE a python/R/Perl script to connect DOOR and download the .opr file

#3.download SRR and TO .fastq
def fastqdump(argv,x):
	# x = 1 means paired-end and x = 0 means single-end
	if x == 1:
		print "running....\n./fastq-dump "+argv+" -split-files -O "+argv+"_input/"
		os.system("./fastq-dump "+argv+" -split-files -O "+argv+"_input/")
	else:
		print "running....\n./fastq-dump "+argv+" -O "+argv+"_input/"
		os.system("./fastq-dump "+argv+" -O "+argv+"_input/")
#4.considering the quality of the FASTQ data

#5.opting a software artifically or automatically

#CHOOSE Rockhopper
def rockhopper(argv,x,y = 0):
	# x = 1 means paired-end and x = 0 means single-end
	# y = 1 means need the SAM file
	os.system("mkdir "+argv+"_output/rockhopper")
	print "running Rockhopper...."
	if y == 0:
		if x == 1:
			os.system("java -Xmx4g -cp Rockhopper.jar Rockhopper -g "+argv+"_ref/ "+argv+"_input/"+argv+"_1.fastq%"+argv+"_input/"+argv+"_2.fastq -o "+argv+"_output/rockhopper -TIME")
		else:
			os.system("java -Xmx4g -cp Rockhopper.jar Rockhopper -g "+argv+"_ref/ "+argv+"_input/"+argv+".fastq -o "+argv+"_output/rockhopper -TIME")
	else:
		if x == 1:
			os.system("java -Xmx4g -cp Rockhopper.jar Rockhopper -g "+argv+"_ref/ "+argv+"_input/"+argv+"_1.fastq%"+argv+"_input/"+argv+"_2.fastq -o "+argv+"_output/rockhopper -TIME -SAM")
			os.system("mv "+argv+"_output/rockhopper/"+argv+"_1.sam "+argv+"_output/"+argv+".sam ")
		else:
			os.system("java -Xmx4g -cp Rockhopper.jar Rockhopper -g "+argv+"_ref/ "+argv+"_input/"+argv+".fastq -o "+argv+"_output/rockhopper -TIME -SAM")
			os.system("mv "+argv+"_output/rockhopper/"+argv+".sam "+argv+"_output/"+argv+".sam ")
	#demonstrating the result
	# os.system("less "+argv+"_output/rockhopper/*operons.txt")

def segemhl(argv,x):
	# output the SAM file
	os.system("./segemhl.x -t 4 -x ref.idx -d "+argv+"_ref/*.fna")
	if x == 1:
		os.system("./segemhl.x -t 4 -i ref.idx -d "+argv+"_ref/*.fna -q "+argv+"_input/"+argv+"_1.fastq -p "+argv+"_input/"+argv+"_2.fastq > "+argv+"_output/"+argv+".sam")
	else:
		os.system("./segemhl.x -t 4 -i ref.idx -d "+argv+"_ref/*.fna -q "+argv+"_input/"+argv+".fastq > "+argv+"_output/"+argv+".sam")
	os.system("rm ref.idx")

def samtools(argv):
	os.system("samtools view -b "+argv+"_output/"+argv+".sam -o "+argv+"_output/"+argv+".bam")
	## reverse
	os.system("samtools view -h -f 16 "+argv+"_output/"+argv+".bam -o "+argv+"_output/"+argv+"_rev.bam")
	## forward
	os.system("samtools view -h -F 16 "+argv+"_output/"+argv+".bam -o "+argv+"_output/"+argv+"_forw.bam")
	## sort
	os.system("samtools sort -o "+argv+"_output/"+argv+"_rev_sort.bam "+argv+"_output/"+argv+"_rev.bam")
	os.system("samtools sort -o "+argv+"_output/"+argv+"_forw_sort.bam "+argv+"_output/"+argv+"_forw.bam")
	## index
	os.system("samtools index "+argv+"_output/"+argv+"_rev_sort.bam")
	os.system("samtools index "+argv+"_output/"+argv+"_forw_sort.bam")
	## count
	os.system("samtools depth -a "+argv+"_output/"+argv+"_forw_sort.bam "+argv+"_output/"+argv+"_rev_sort.bam > "+argv+"_output/"+argv+"_count")
#CHOOSE CONDOP
def CONDOP(argv,x,align = "segemhl"):
	#align
	if align == "segemhl":
		segemhl(argv,x)
	else:
		rockhopper(argv,x,y=1)
	#samtools:conut coverage depth
	samtools(argv)
	os.system("python extract_col_3_4.py -i "+argv+"_count -o "+argv+"_table")
	os.system("Rscript CONDOP_script.R "+argv)
	os.system("mv COP.CONDOP.txt "+argv+"_output/CONDOP/"+argv+"_operons.txt")
	#demonstrating the result
	os.system("less "+argv+"_output/CONDOP/*operons.txt")

#CHOOSE RNAseg
def RNAseg(argv,x):
	segemhl(argv,x)
	if x == 0:
		os.system("python sam2grp.py "+argv+"_output/"+argv+".sam")
		#os.system("RNAseg -t 8 -f "+argv+"fwd,grp")
	if x == 1:
		insert_size = ""
		print "cannot WRITE the code NOW!\ncannot know the insert_size"
		# os.system("python sam2grp.py "+argv+"_output/"+argv+".sam "+insert_size)
		#os.system("RNAseg -t 8 -f "+argv+"fwd,grp")

if __name__ == '__main__':
	if len(sys.argv)==1:
		usage()
	else:
		arg = sys.argv[1]
		checkargv(arg)
		# print "OK"
		information={'Organism':"",'Instrument':"",'Layout':""}
		getinformation(arg,information)		
		x = paired_or_single(information["Layout"])
		print x
		# x = 1 means paired-end and x = 0 means single-end
		mkdir(arg)
		fastqdump(arg,x)
		download_annotion(arg,information['Organism'])
		rockhopper(arg,x)
		#CONDOP(arg,x)
		#RNAseg(arg,x)
