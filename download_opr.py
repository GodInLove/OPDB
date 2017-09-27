# -*- coding: UTF-8 -*-
import sys
import getopt
import urllib2
import re
import urllib

def usage():
	print "Please use python download_opr.py organism_name ref_path"

def Schedule(a,b,c):
	per = 100.0 * a * b / c
	if per > 100 :
		per = 100
	print '%.2f%%' % per

def search_opr(organism_name):
	root_url = "http://csbl.bmb.uga.edu/DOOR/"
	# browser search string
	pat_init = re.compile(r"\s")
	serach_str = re.sub(pat_init,"%20",organism_name)
	# try connect DOOR
	retry = -10
	while retry < 0:	
		try:
			url = "http://csbl.bmb.uga.edu/DOOR/search_ajax.php?keyword="+serach_str+"&mode=DataTable&sEcho=1&iColumns=6&sColumns=&iDisplayStart=0&iDisplayLength=15&mDataProp_0=0&mDataProp_1=1&mDataProp_2=2&mDataProp_3=3&mDataProp_4=4&mDataProp_5=5&_=1506424463060"
			print "search in the DOOR database,please wait...\n"
			page = urllib2.urlopen(url,timeout=10)
			html = page.read()
			#check whether the opr file is available or not
			pat_notfind = re.compile(r'iTotalRecords\"\:([0-9]+)')
			searchobj = re.findall(pat_notfind,html)
			if searchobj[0] == '0':
				print "the DOOR database has not relative .opr file,please change another method"
				sys.exit(2)
			# when the opr file is available
			_str = html.split("[")
			pat1 = re.compile(r'(operon\.php\?id\=[0-9]+).*'+organism_name.split()[-1])
			for item in _str:
				result1 = re.findall(pat1,item)
				if result1:
					break
			# print result1
		except:
			retry+=1
			if retry == 0:
				print "wrong!!Cannot connect DOOR"
				sys.exit(2)
	#the next page	
	retry = -10
	while retry < 0:
		try:	
			url1 = root_url + result1[0]
			print "search in the next page,please wait...\n"
			page = urllib2.urlopen(url1,timeout=20)
			html = page.read()
			pat2 = re.compile(r'\<li\>\<a\shref\=\"(.*)\".*\<\/li\>\s+\<li\>\<a\shref\=\"operon\.php\?id\=3492')
			result2 = re.findall(pat2,html)
			if not result2:
				print "cannot find the url,please contact with the author"
			# print result2
		except:
			retry+=1
			if retry == 0:
				print "wrong!!Cannot connect DOOR"
				sys.exit(2)
	#the next page
	retry = -10
	while retry < 0:
		try:
			url2 = root_url + result2[0]
			page = urllib2.urlopen(url2,timeout=20)
			html = page.read()
			pat3 = re.compile(r'window\.location\=\'(downloadNCoperon.php.*)\'\"\svalue\=\"Download\soperon\stable\"\>')
			result3 = re.findall(pat3,html)
			if not result2:
				print "cannot find the url,please contact with the author"
		except:
			retry+=1
			if retry == 0:
				print "wrong!!Cannot connect DOOR"
				sys.exit(2)
	#return the finall url
	url_download = root_url + result3[0]
	return url_download

def download_opr(url_download,ref_path)
	pat_name = re.compile(r'[0-9]+')
	result = re.findall(pat_name,url_download)
	local = ref_path +"/"+ result[0] + ".opr"
	urllib.urlretrieve(url_download,local,Schedule)


if __name__ == '__main__':
	if len(sys.argv) < 1:
		usage()
		sys.exit()
	organism_name = sys.argv[1]
	ref_path = sys.argv[2]
	url = search_opr(organism_name)
	print "the download link is:\n"+url
	download_opr(url,ref_path)
