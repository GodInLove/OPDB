import os



def srr_n_to_fastq(srr_n, layout, output_path):
    output_path = output_path + srr_n
    if not os.path.exists(output_path):
        os.system("mkdir " + output_path)
    if layout == 1:
        print("running....\nfastq-dump " + srr_n + " -split-files -O " + output_path)
        if not os.path.exists(output_path + "/" +srr_n + "_1.fastq"):
            os.system("tools/fastq-dump " + srr_n + " -split-files -O " + output_path)
    else:
        print("running....\nfastq-dump " + srr_n + " -O " + output_path)
        if not os.path.exists(output_path + "/" + srr_n + ".fastq"):
            os.system("tools/fastq-dump " + srr_n + " -O " + output_path)
