import os


# _dir is a list with (input_path, ref_path, output_path)
def segemehl_align(srr_n, x, _dir, process_n):
    # output the SAM file
    os.system("tools/segemehl.x -t " + process_n + " -x ref.idx -d " + _dir[1] + "/*.fna")
    if x == 1:
        os.system("tools/segemehl.x -t " + process_n + " -i ref.idx -d " + _dir[1] + "/*.fna -q " + _dir[
            0] + "/" + srr_n + "_1.fastq -p " + _dir[0] + "/" + srr_n + "_2.fastq > " + _dir[2] + "/" + srr_n + ".sam")
    else:
        os.system("tools/segemehl.x -t " + process_n + " -i ref.idx -d " + _dir[1] + "/*.fna -q " + _dir[
            0] + "/" + srr_n + ".fastq > " + _dir[2] + "/" + srr_n + ".sam")
    os.system("rm ref.idx")
