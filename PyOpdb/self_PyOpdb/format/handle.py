import re

import os


def extract_Synonym_sub(infile, ref):
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


def extract_Synonym(_dir, srr_n):
    synonym = []
    for item in os.listdir(_dir[1]):
        if ".ptt" in item:
            ptt_file = item
            synonym = extract_Synonym_sub(_dir[1] + "/" + ptt_file, synonym)
        if ".rnt" in item:
            rnt_file = item
            synonym = extract_Synonym_sub(_dir[1] + "/" + rnt_file, synonym)
    # print(len(synoym))
    f = open(_dir[2] + "/_operons.txt", 'r')
    content = f.read()
    f.close()
    for it in synonym:
        content = content.replace(it[0], it[1])
        f = open(_dir[2] + "/" + "_operon.txt", 'w')
        f.write(content)
        f.close()
    os.system("rm " + _dir[2] + "/_operons.txt")


def extract_wig(wig_path, gff_path):
    return 0


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
