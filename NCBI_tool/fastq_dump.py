import os


def sra2fastq(srr_n, x, input_path):
    # x = 1 means paired-end and x = 0 means single-end
    if x == 1:
        print("running....\nfastq-dump " + srr_n + " -split-files -O " + input_path)
        os.system("tools/fastq-dump " + srr_n + " -split-files -O " + input_path)
    else:
        print("running....\nfastq-dump " + srr_n + " -O " + input_path)
        os.system("tools/fastq-dump " + srr_n + " -O " + input_path)
