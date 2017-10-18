import os

from srr2sam_tool.rockhopper_align import rockhopper_align
from srr2sam_tool.segemehl import segemehl_align
from format_handle.extract_col import extract_3_4, extract_some
from format_handle.samtool_handle import samtools


def CONDOP_operon_predict(srr_n, x, _dir, process_n, align="segemehl"):
    # align
    if align == "segemehl":
        segemehl_align(srr_n, x, _dir, process_n)
    else:
        rockhopper_align(srr_n, x, _dir, process_n)
    # samtools:conut coverage depth
    samtools(srr_n, _dir[2])
    extract_3_4(_dir[2] + "/" + srr_n + "_count", _dir[2] + "/" + srr_n + "_table")
    os.system("rm " + _dir[2] + "/" + srr_n + "_count")
    os.system("mkdir " + _dir[2] + "/CONDOP")
    os.system("Rscript scripts/bin/CONDOP_script.R " + srr_n + " " + _dir[1] + " " + _dir[2])
    extract_some(_dir[2] + "/CONDOP/" + "COP.CONDOP.1.txt",_dir[2] + "/CONDOP/" + srr_n + "_operons.txt")
    print("Operons written to file:\t" + _dir[2] + "/CONDOP/" + srr_n + "_operon.txt")
    # os.system("rm "+_dir[2] + "/COP.CONDOP.txt")
