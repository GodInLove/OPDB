import os


# _dir is a list with (input_path, ref_path, output_path)
def rockhopper_align(srr_n, x, _dir, process_n):
    if x == 1:
        os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
            0] + "/" + srr_n + "_1.fastq%" + _dir[0] + "/" + srr_n + "_2.fastq -o " + _dir[
                      2] + "/rockhopper -TIME -SAM")
        os.system("mv " + _dir[2] + "/rockhopper/" + srr_n + "_1.sam " + _dir[2] + "/" + srr_n + ".sam")
    else:
        os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
            0] + "/" + srr_n + ".fastq -o " + _dir[2] + "/rockhopper -TIME -SAM")
        os.system("mv " + _dir[2] + "/rockhopper/" + srr_n + ".sam " + _dir[2] + "/" + srr_n + ".sam")
