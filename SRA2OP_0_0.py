#!bin/python2.7
import os
import getopt
import sys
import re
import urllib2

# Basic V0.0

def usage():
	print "\nUseage:python SRA2OP_0_0.py -i SRR5486953"
	print "\nOr python SRA2OP_0_0.py -i SRR_number -o output_path -m <int> -p <int> "
	print "\n-i|--input[string]\tSRRnumber(default:null)"
	print "\n-o|--output[string]\toutput_dir(default:current dir)"
	print "\n-m|--method[int]\tchoose the software, (0,1,2) means (rockhopper,CONDOP,RNAseg)(default:0)"
	print "\n-p|--processor[int]\tset the number of processor(default:4)"

#1.input SRRXXXX and check it
def checksrr_n(srr_n):
	pat = re.compile(r'[SRA]+[0-9]+')
	searchobj = re.search(pat,srr_n)
	if not searchobj:
		print "wrong SRAnumber!"
		usage()
		sys.exit(2)

def check_SRR(information):
	#check whether the SRR data is up to standard or not
	if information['Organism'] == "Escherichia coli":
		not_qualified()
	else:
		pass

def checkprocess_n(process_n):
	if int(process_n) < 0:
		print "wrong process number!"
		usage()
		sys.exit(2)
	elif int(process_n) >100:
		print "wrong process number!"
		usage()
		sys.exit(2)

def checkmethod(method):
	if int(method) < 0:
		print "wrong method number!"
		usage()
		sys.exit(2)
	elif int(method) > 2:
		print "wrong method number!"
		usage()
		sys.exit(2)

def checkoutput_path(output_path):
	if len(output_path) > 0:
		if output_path[-1] != "/":
			output_path = output_path+"/"
	return output_path

def mkdir(srr_n,output_path):
	input_dir = output_path+srr_n+"_input"
	ref_dir = output_path+srr_n+"_ref"
	output_dir = output_path+srr_n+"_output"
	os.system("mkdir "+output_path+srr_n+"_input")
	os.system("mkdir "+output_path+srr_n+"_ref")
	os.system("mkdir "+output_path+srr_n+"_output")
	return [input_dir,ref_dir,output_dir]

def not_qualified():
	print "the SRR data is not qualified!"
	sys.exit(2)

#2.figure out something and prepare reference data 
'''
figuring out different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
figuring out the parameter settings of the software
figuring out whether the data is paired-end or not
figuring out .opr file
figuring out .gff .rnt .fna .gtf
'''
#WRITE a python/R/Perl script to connect NCBI API and get those information
def getinformation(srr_n,information):
	retry = -10
	while retry < 0:
		# print retry
		try:
			url = 'https://www.ncbi.nlm.nih.gov/sra/?term='+srr_n+'%5BAll+Fields%5D+AND+"biomol+rna"%5BProperties%5D'
			print "open link:\n"+url+"\nplease wait...\n"
			page = urllib2.urlopen(url,timeout=10)
			html = page.read()
			retry = 1
			pat_notfind = re.compile(r'No items found')
			searchobj = re.search(pat_notfind,html)
			if searchobj:
				print "SRA Experiment "+srr_n+" is not public or not RNA. please change your SRAnumber"
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
			for item in information:
				print(item+": "+information[item]+"\n")
			# f = open(srr_n+"_information.txt",'w')
			# print "\nNow you can open your workspace, and you will find a file named "+srr_n+"_information.txt. It stores some BASIC information about your SRR data."
			# for item in information:
			# 	f.write(item+": "+information[item]+"\n")
			# f.close()
		except:
			retry+=1
			if retry == 0:
				print "wrong!!Cannot connect NCBI"
				sys.exit(2)

#WRITE a python/R/Perl script to connect NCBI API and download the ref.fna ref.gff ref.rnt ref.gtf
def download_annotion(srr_n,organism,_dir):
	print "downloading annotion:.gtf .rnt .gff......."
	# assume that finished it
	# download the annotion in the ref_dir
	os.system("cp NC_000913_ref/* "+_dir[1])
	print organism+"annotion files were downloaded in the "+_dir[1]

#WRITE a python/R/Perl script to figure out whether it is paired-end or not
def paired_or_single(layout):
	pat_paired = re.compile(r'[padPAD]')
	searchobj = re.search(pat_paired,layout)
	if searchobj:
		return 1
	else:
		return 0



#3.download SRR and TO .fastq
def fastqdump(srr_n,x,_dir):
	# x = 1 means paired-end and x = 0 means single-end
	if x == 1:
		print "running....\n./fastq-dump "+srr_n+" -split-files -O "+_dir[0]
		os.system("./fastq-dump "+srr_n+" -split-files -O "+_dir[0])
	else:
		print "running....\n./fastq-dump "+srr_n+" -O "+_dir[0]
		os.system("./fastq-dump "+srr_n+" -O "+_dir[0])
#4.considering the quality of the FASTQ data

#5.opting a software artifically or automatically

#CHOOSE Rockhopper
def rockhopper(srr_n,x,_dir,process_n,y = 0):
	# x = 1 means paired-end and x = 0 means single-end
	# y = 1 means need the SAM file
	os.system("mkdir "+_dir[2]+"/rockhopper")
	print "running Rockhopper...."
	if y == 0:
		if x == 1:
			os.system("java -Xmx3g -cp Rockhopper.jar Rockhopper -p "+process_n+" -g "+_dir[1]+" "+_dir[0]+"/"+srr_n+"_1.fastq%"+_dir[0]+"/"+srr_n+"_2.fastq -o "+_dir[2]+"/rockhopper -TIME")
		else:
			os.system("java -Xmx3g -cp Rockhopper.jar Rockhopper -p "+process_n+" -g "+_dir[1]+" "+_dir[0]+"/"+srr_n+".fastq -o "+_dir[2]+"/rockhopper -TIME")
	else:
		if x == 1:
			os.system("java -Xmx3g -cp Rockhopper.jar Rockhopper -p "+process_n+" -g "+_dir[1]+" "+_dir[0]+"/"+srr_n+"_1.fastq%"+_dir[0]+"/"+srr_n+"_2.fastq -o "+_dir[2]+"/rockhopper -TIME -SAM")
			os.system("mv "+_dir[2]+"/rockhopper/"+srr_n+"_1.sam "+_dir[2]+"/"+srr_n+".sam")
		else:
			os.system("java -Xmx3g -cp Rockhopper.jar Rockhopper -p "+process_n+" -g "+_dir[1]+" "+_dir[0]+"/"+srr_n+".fastq -o "+_dir[2]+"/rockhopper -TIME -SAM")
			os.system("mv "+_dir[2]+"/rockhopper/"+srr_n+".sam "+_dir[2]+"/"+srr_n+".sam")
	#demonstrating the result
	# os.system("less "+srr_n+"_output/rockhopper/*operons.txt")

def segemehl(srr_n,x,_dir,process_n):
	# output the SAM file
	os.system("./segemehl.x -t "+process_n+" -x ref.idx -d "+_dir[1]+"/*.fna")
	if x == 1:
		os.system("./segemehl.x -t "+process_n+" -i ref.idx -d "+_dir[1]+"/*.fna -q "+_dir[0]+"/"+srr_n+"_1.fastq -p "+_dir[0]+"/"+srr_n+"_2.fastq > "+_dir[2]+"/"+srr_n+".sam")
	else:
		os.system("./segemehl.x -t "+process_n+" -i ref.idx -d "+_dir[1]+"/*.fna -q "+_dir[0]+"/"+srr_n+".fastq > "+_dir[2]+"/"+srr_n+".sam")
	os.system("rm ref.idx")

def samtools(srr_n,_dir):
	os.system("samtools view -b "+_dir[2]+"/"+srr_n+".sam -o "+_dir[2]+"/"+srr_n+".bam")
	## reverse
	os.system("samtools view -h -f 16 "+_dir[2]+"/"+srr_n+".bam -o "+_dir[2]+"/"+srr_n+"_rev.bam")
	## forward
	os.system("samtools view -h -F 16 "+_dir[2]+"/"+srr_n+".bam -o "+_dir[2]+"/"+srr_n+"_forw.bam")
	## sort
	os.system("samtools sort -o "+_dir[2]+"/"+srr_n+"_rev_sort.bam "+_dir[2]+"/"+srr_n+"_rev.bam")
	os.system("samtools sort -o "+_dir[2]+"/"+srr_n+"_forw_sort.bam "+_dir[2]+"/"+srr_n+"_forw.bam")
	## index
	os.system("samtools index "+_dir[2]+"/"+srr_n+"_rev_sort.bam")
	os.system("samtools index "+_dir[2]+"/"+srr_n+"_forw_sort.bam")
	## count
	os.system("samtools depth -a "+_dir[2]+"/"+srr_n+"_forw_sort.bam "+_dir[2]+"/"+srr_n+"_rev_sort.bam > "+_dir[2]+"/"+srr_n+"_count")
#CHOOSE CONDOP
def CONDOP(srr_n,x,_dir,process_n,align = "segemehl"):
	#align
	if align == "segemehl":
		segemehl(srr_n,x,_dir,process_n)
	else:
		rockhopper(srr_n,x,_dir,process_n,y=1)
	#samtools:conut coverage depth
	samtools(srr_n,_dir)
	os.system("python extract_col_3_4.py -i "+_dir[2]+"/"+srr_n+"_count -o "+_dir[2]+"/"+srr_n+"_table")
	os.system("rm "+_dir[2]+"/"+srr_n+"_count")
	os.system("Rscript CONDOP_script.R "+srr_n+" "+_dir[1]+" "+_dir[2])
	os.system("mv "+_dir[2]+"/COP.CONDOP.txt "+_dir[2]+"/CONDOP/"+srr_n+"_operons.txt")
	#demonstrating the result
	os.system("less "+_dir[2]+"/CONDOP/*operons.txt")

#CHOOSE RNAseg
def RNAseg(srr_n,x,_dir,process_n):
	segemehl(srr_n,x,_dir,process_n)
	if x == 0:
		os.system("python sam2grp.py "+_dir[2]+"/"+srr_n+".sam")
		#os.system("RNAseg -t 8 -f "+srr_n+"fwd,grp")
	if x == 1:
		insert_size = ""
		print "cannot WRITE the code NOW!\ncannot know the insert_size"
		# os.system("python sam2grp.py "+srr_n+"_output/"+srr_n+".sam "+insert_size)
		#os.system("RNAseg -t 4 -f "+srr_n+"fwd,grp")
def main(argv):
	srr_n = ""
	method = 0
	process_n = 4
	output_path = ""
	try:
		opts, args = getopt.getopt(argv[1:],"hi:o:m:p:",["help","input=","output=","method=","processor="])
		for opt,arg in opts:
			if opt in ("-h","--help"):
				usage()
				sys.exit()
			elif opt in ("-i","--input"):
				srr_n = arg
			elif opt in ("-o","--output"):
				output_path = arg
			elif opt in ("-m","--method"):
				method = int(arg)
			elif opt in ("-p","--processor"):
				process_n = int(arg)
		checksrr_n(srr_n)
		checkmethod(method)
		checkprocess_n(process_n)
		output_path = checkoutput_path(output_path)
		print "OK"
		information={'Organism':"",'Instrument':"",'Layout':""}
		getinformation(srr_n,information)		
		x = paired_or_single(information["Layout"])
		print x
		# x = 1 means paired-end and x = 0 means single-end
		_dir = mkdir(srr_n,output_path)
		#_dir has 3 dir, input_dir, ref_dir, out_dir
		fastqdump(srr_n,x,_dir)
		download_annotion(srr_n,information['Organism'],_dir)
		if method == 0:
			rockhopper(srr_n,x,_dir,str(process_n))
		elif method == 1:
			#WRITE a python/R/Perl script to connect DOOR and download the .opr file
			os.system("python download_opr.py "+information['Organism']+" "+_dir[1])
			CONDOP(srr_n,x,_dir,str(process_n))
		elif method == 2:
			pass
			#RNAseg(srr_n,x,_dir,str(process_n))
	except getopt.GetoptError:
		usage()
		sys.exit(2)

if __name__ == '__main__':
	main(sys.argv)