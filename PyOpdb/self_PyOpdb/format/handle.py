import re

import os

import sys

import NCBI_API


def extract_Synonym(ref_path, result_path):
    """
    :method: a method to change like "b0001" to "thrL"
    :param ref_path: string
    :param result_path: string
    :return:void(some files in result_path)
    """
    synonym = []
    for item in os.listdir(ref_path):
        if ".ptt" in item:
            ptt_file = item
            synonym = extract_Synonym_sub(ref_path + "/" + ptt_file, synonym)
        if ".rnt" in item:
            rnt_file = item
            synonym = extract_Synonym_sub(ref_path + "/" + rnt_file, synonym)
    # print(len(synoym))
    for item in os.listdir(result_path):
        if "operon" in item:
            operon_file = item
    f = open(result_path + "/" + operon_file, 'r')
    content = f.read()
    f.close()
    # change like "b0001" to "thrL"
    for it in synonym:
        content = content.replace(it[0], it[1])
        f = open(result_path + "/" + "_operon.txt", 'w')
        f.write(content)
        f.close()
    os.system("rm " + result_path + "/*_operons.txt")
    print("Operons written to file:\t" + result_path + "_operon.txt")


def extract_wig(wig_path, gff_path, result_path):
    """
    :method: a method to convert wig to bigwig
    :param wig_path: string
    :param gff_path: string
    :param result_path: string
    :return: void(somefiles in result_path)
    """
    # get chrome.sizes
    chrome_path = result_path + "/chrome.size"
    f = open(gff_path, 'r')
    content = f.read()
    f.close()
    chrome = NCBI_API.findall_pat("\#\#sequence\-region ([NC\_0-9\.]+ 1 [0-9]+)", content)
    chrome = chrome[0].split(" ")
    f = open(chrome_path, 'w')
    f.write(chrome[0] + "\t" + chrome[2] + "\n")
    f.close()
    # remove the "track_name" in wig file
    f = open(wig_path, 'r')
    trash = f.readline()
    content = f.read()
    f.close()
    os.system("rm " + wig_path)
    f = open(wig_path, 'w')
    f.write(content)
    f.close()
    # run the bash to convert wig to bigwig
    tools_path = "/home/yaodongliu/OPDB/PyOpdb/tools/wigToBigwig"
    # print(tools_path + " " + wig_path + " " + chrome_path + " " + result_path + "/operon.bw")
    os.popen(tools_path + " " + wig_path + " " + chrome_path + " " + result_path + "/operon.bw")


def extract_Synonym_sub(infile, ref):
    """
    :method: a sub-method
    :param infile: string(path)
    :param ref: string
    :return:
    """
    f = open(infile, 'r')
    for i in range(0, 3):
        f.readline()
    while True:
        line = f.readline().strip()
        if not line:
            break
        li = line.split("\t")
        ref.append([li[4], li[5]])
    f.close()
    return ref


def extract_3_4(infile, outfile):
    f = open(infile, 'r')
    out_f = open(outfile, 'w')
    out_f.write("fwd,rev\n")
    while True:
        line = f.readline().strip()
        if not line:
            break
        li = line.split()
        out_f.write(li[2] + "," + li[3] + "\n")

    f.close()
    out_f.close()


def extract_some(infile, outfile):
    f = open(infile, 'r')
    out_f = open(outfile, 'w')
    out_f.write("Start\tStop\tStrand\tNumber of Genes\tGenes\n")
    line = f.readline()
    pat = re.compile(r'\-')
    while True:
        line = f.readline().strip()
        if not line:
            break
        li = line.split("\t")
        li[7] = re.sub(pat, ",", li[7])
        out_f.write(li[2] + "\t" + li[3] + "\t" + li[1] + "\t" + li[4] + "\t" + li[7] + "\n")

    f.close()
    out_f.close()
