import urllib.request

import re

import sys


def getinformation(srr_n, information):
    try:
        url = 'https://www.ncbi.nlm.nih.gov/sra/?term=' + srr_n + '%5BAll+Fields%5D+AND+"biomol+rna"%5BProperties%5D'
        print("search in the SRA database,please wait...\n")
        req = urllib.request.Request(url)
        response = urllib.request.urlopen(req)
        html = response.read()
        html = html.decode('utf-8')
        pat_nothing = re.compile(r'No items found')
        searchable = re.search(pat_nothing, html)
        if searchable:
            print("SRA Experiment " + srr_n + " is not public or not RNA. please change your SRAnumber")
            sys.exit(2)
        # find Organism
        pat_org = re.compile(r'Organism\: \<span\>\<a href\=.*\"\>([a-zA-Z0-9\.\s\-]+).*expand showed sra-full-data')
        org = re.findall(pat_org, html)
        if len(org) == 1:
            information['Organism'] = org[0]
        else:
            print("Cannot find the Organism.WHAT happen?the result is", org)
            sys.exit(2)
        # find Instrument
        pat_ins = re.compile(r'Instrument\: \<span\>([A-Za-z0-9\s]+)\<\/span\>\<\/div\>\<div\>Strategy\:')
        ins = re.findall(pat_ins, html)
        if len(ins) == 1:
            information['Instrument'] = ins[0]
        else:
            print("Cannot find the Instrument.WHAT happen?the result is", ins)
            sys.exit(2)
        # find Layout
        pat_lay = re.compile(r'Layout: <span>([PAIREDSINGLEpairedsinglenN]+).*sra-full-data')
        lay = re.findall(pat_lay, html)
        if len(lay) == 1:
            information['Layout'] = lay[0]
        else:
            print("Cannot find the Layout.WHAT happen?the result is", org)
            sys.exit(2)
        for item in information:
            print(item + ": " + information[item] + "\n")
        # f = open(srr_n+"_information.txt",'w')
        # print "\nNow you can open your workspace, and you will find a file named "+srr_n+"_information.txt. It stores some BASIC information about your SRR data."
        # for item in information:
        # 	f.write(item+": "+information[item]+"\n")
        # f.close()
    except:
        print("wrong!!Cannot connect NCBI")
        sys.exit(2)


# WRITE a python/R/Perl script to figure out whether it is paired-end or not
def paired_or_single(layout):
    pat_paired = re.compile(r'[padPAD]')
    searchable = re.search(pat_paired, layout)
    if searchable:
        return 1
    else:
        return 0
