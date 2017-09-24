# -*- coding: UTF-8 -*-
import sys
import getopt

def usage():
	print "Please use python extract_col_3_4.py -i input -o output"

def extract(infile,outfile):
	f=open(infile,'r')
	out_f = open(outfile,'w')
	out_f.write("fwd,rev\n")
	while True:
		line = f.readline().strip()
		if not line:
			break
		li = line.split()
		out_f.write(li[2]+","+li[3]+"\n")

	f.close()
	out_f.close()

def main(argv):
	inputfile = ""
	outputfile = ""
	try:
		opts, args = getopt.getopt(argv[1:],"hi:o:",["help","input=","output="])
		for opt,arg in opts:
			if opt in ("-h","--help"):
				usage()
				sys.exit()
			elif opt in ("-i","--input"):
				inputfile = arg
			elif opt in ("-o","--output"):
				outputfile = arg
		extract(inputfile,outputfile)
	except getopt.GetoptError:
		usage()
		sys.exit(2)

if __name__ == "__main__":
   main(sys.argv)