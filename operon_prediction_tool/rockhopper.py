import os


# _dir is a list with (input_path, ref_path, output_path)
import sys

from format_handle.extract_col import extract_Synonym
from format_handle.samtool_handle import samtools


def rockhopper_operon_predict(srr_n, x, _dir, process_n):
    if os.path.exists(_dir[1] + "/*.gff"):
        os.system("rm " + _dir[1] + "/*.gff")
    # x = 1 means paired-end and x = 0 means single-end
    if not os.path.exists(_dir[2]+"/rockhopper"):
        os.system("mkdir " + _dir[2] + "/rockhopper")
        print("running Rockhopper....")
        if x == 1:
            os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
                0] + "/" + srr_n + "_1.fastq%" + _dir[0] + "/" + srr_n + "_2.fastq -o " + _dir[
                          2] + "/rockhopper -TIME -SAM")
        else:
            os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
                0] + "/" + srr_n + ".fastq -o " + _dir[2] + "/rockhopper -TIME -SAM")
        extract_Synonym(_dir, srr_n)
        samtools(srr_n, _dir[2] + "/rockhopper")
        os.system("rm " + _dir[2] + "/rockhopper/genomeBrowserFiles/_diff*")
        os.system("mv " + _dir[2] + "/rockhopper/genomeBrowserFiles/*wig " + _dir[2])
        os.system("mv " + _dir[2] + "/rockhopper/*bw " + _dir[2])
        print("Operons written to file:\t" + _dir[2] + "/rockhopper/" + srr_n + "_operon.txt")
    else:
        print("\nthe result has done ! Please check the path： "+_dir[2]+"/rockhopper\n")

