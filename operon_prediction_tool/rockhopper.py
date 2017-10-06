import os


# _dir is a list with (input_path, ref_path, output_path)
import sys

from format_handle.extract_col import extract_Synonym


def rockhopper_operon_predict(srr_n, x, _dir, process_n):
    os.system("rm " + _dir[1] + "/*.gff")
    # x = 1 means paired-end and x = 0 means single-end
    if not os.path.exists(_dir[2]+"/rockhopper"):
        os.system("mkdir " + _dir[2] + "/rockhopper")
    else:
        print("\nthe result has done ! Please check the path： "+_dir[2]+"/rockhopper\n")
        sys.exit(2)
    print("running Rockhopper....")
    if x == 1:
        os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
            0] + "/" + srr_n + "_1.fastq%" + _dir[0] + "/" + srr_n + "_2.fastq -o " + _dir[2] + "/rockhopper -TIME")
    else:
        os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
            0] + "/" + srr_n + ".fastq -o " + _dir[2] + "/rockhopper -TIME")
    extract_Synonym(_dir,srr_n)
    print("Operons written to file:\t" + _dir[2] + "/rockhopper/" + srr_n + "_operon.txt")
