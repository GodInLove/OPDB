#!bin/python2.7
import os
import getopt
import sys
import re

# Basic V0.0

def usage():
	print "Useage:python SRA2OP_0_0.py SRR0000"

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

#2.figure out something and prepare reference data 
'''
figuring out different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
figuring out the parameter settings of the software
figuring out whether the data is paired-end or not
figuring out .opr file
figuring out .gff .rnt .fna .gtf
'''
#WRITE a python/R/Perl script to connect NCBI API and get those information

#WRITE a python/R/Perl script to connect NCBI API and download the ref.fna ref.gff ref.rnt ref.gtf

#WRITE a python/R/Perl script to connect DOOR and download the .opr file

#3.download SRR and TO .fastq
def fastqdump(argv,x):
	# x 1 means paired-end, 0 means single-end
	if x==1:
		os.system("./fastq-dump "+argv+" -O "+argv+"_input/")
	elif x==0:
		os.system("./fastq-dump $SRRnumber -split-files -O "+argv+"_input/")
	else:
		pass

#4.considering the quality of the FASTQ data

#5.opting a software artifically or automatically

#example with Rockhopper
def rockhopper(argv):
	os.system("mkdir "+argv+"_output/rockhopper")
	#already
	os.system("java -Xmx4g -cp Rockhopper.jar Rockhopper -g "+argv+"_ref/ "+argv+"_input/"+argv+".fastq -o "+argv+"_output/rockhopper -TIME")
	#os.system("java -Xmx4g -cp Rockhopper.jar Rockhopper -g "+argv+"_ref/ "+argv+"_input/"+argv+"_1.fastq%"+argv+"_input/"+argv+"_2.fastq -o Rockhopper_"+argv+"_output/ -TIME")
	#demonstrating the result
	os.system("cat "+argv+"_output/rockhopper/*operons.txt")

if __name__ == '__main__':
	if len(sys.argv)==1:
		usage()
	else:
		arg = sys.argv[1]
		checkargv(arg)
		# print "OK"