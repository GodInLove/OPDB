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


def extract_Synonym(_dir, srr_n):
    for item in os.listdir(_dir[1]):
        if ".ptt" in item:
            ptt_file = item
        if ".rnt" in item:
            rnt_file = item
    ref = []
    ref = extract_Synonym_sub(_dir[1] + "/" + ptt_file, ref)
    ref = extract_Synonym_sub(_dir[1] + "/" + rnt_file, ref)
    # print(len(ref))
    for item in os.listdir(_dir[2] + "/rockhopper"):
        if "operons.txt" in item:
            operon_file = item
            f = open(_dir[2] + "/rockhopper/" + operon_file, 'r')
            content = f.read()
            f.close()
            for it in ref:
                content = content.replace(it[0], it[1])
            f = open(_dir[2] + "/rockhopper/" + srr_n + "_operon.txt", 'w')
            f.write(content)
            f.close()
            os.system("rm " + _dir[2] + "/rockhopper/" + operon_file)
        else:
            pass
