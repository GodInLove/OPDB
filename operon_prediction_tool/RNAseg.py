import os

from align_tool.segemehl import segemehl_align


def RNAseg_operon_predict(srr_n, x, _dir, process_n):
    segemehl_align(srr_n, x, _dir, process_n)
    if x == 0:
        os.system("python scripts/sam2grp.py " + _dir[2] + "/" + srr_n + ".sam")
    # os.system("RNAseg -t 8 -f "+srr_n+"fwd,grp")
    if x == 1:
        # insert_size = ""
        print("cannot WRITE the code NOW!\ncannot know the insert_size")
        # os.system("python sam2grp.py "+srr_n+"_output/"+srr_n+".sam "+insert_size)
        # os.system("RNAseg -t "+process_n+" -f "+srr_n+"fwd,grp")
