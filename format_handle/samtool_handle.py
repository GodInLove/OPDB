import os


def samtools(srr_n,output_path):
    os.system("samtools view -b " + output_path + "/" + srr_n + ".sam -o " + output_path + "/" + srr_n + ".bam")
    # reverse
    os.system("samtools view -h -f 16 " + output_path + "/" + srr_n + ".bam -o " + output_path + "/" + srr_n + "_rev.bam")
    # forward
    os.system("samtools view -h -F 16 " + output_path + "/" + srr_n + ".bam -o " + output_path + "/" + srr_n + "_forw.bam")
    # sort
    os.system("samtools sort -o " + output_path + "/" + srr_n + "_rev_sort.bam " + output_path + "/" + srr_n + "_rev.bam")
    os.system("samtools sort -o " + output_path + "/" + srr_n + "_forw_sort.bam " + output_path + "/" + srr_n + "_forw.bam")
    # index
    os.system("samtools index " + output_path + "/" + srr_n + "_rev_sort.bam")
    os.system("samtools index " + output_path + "/" + srr_n + "_forw_sort.bam")
    # count
    # os.system(
    #     "samtools depth -a " + output_path[2] + "/" + srr_n + "_forw_sort.bam " + output_path[2] + "/" + srr_n + "_rev_sort.bam > " +
    #     output_path + "/" + srr_n + "_count")
    # to bigwig
    os.system("bamCoverage -b " + output_path + "/" + srr_n + "_rev_sort.bam -o " + output_path + "/" + srr_n + "_rev_sort.bw")
    os.system(
        "bamCoverage -b " + output_path + "/" + srr_n + "_forw_sort.bam -o " + output_path + "/" + srr_n + "_forw_sort.bw")