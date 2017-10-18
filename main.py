import getopt
import os
import sys

from NCBI_tool.getinformation import getinformation, paired_or_single

from check_input.check_argv import usage, check_args
from srr2sam_tool.fastq_dump import sra2fastq
from download_annotion.download_annot import download_annotion
from download_annotion.download_opr import download_opr
from format_handle.result2jbrowse import res2jbrowse
from operon_prediction_tool.CONDOP import CONDOP_operon_predict
from operon_prediction_tool.rockhopper import rockhopper_operon_predict


def makedir(srr_n, output_path, kegg_id):
    input_dir = output_path + srr_n + "_input"
    ref_dir = output_path + kegg_id
    output_dir = output_path + srr_n + "_output"
    if not os.path.exists(input_dir):
        switch_1 = 0
        os.system("mkdir " + input_dir)
    else:
        switch_1 = 1
        print("\n" + input_dir + "\tthe dir exists.\n")
    if not os.path.exists(ref_dir):
        switch_2 = 0
        os.system("mkdir " + ref_dir)
    else:
        switch_2 = 1
        print("\n" + ref_dir + "\tthe dir exists.\n")
    if not os.path.exists("/home/lyd/webapps/JBrowse/"+srr_n):
        switch_3 = 0
        if not os.path.exists(output_dir):
            os.system("mkdir " + output_dir)
    else:
        switch_3 = 1
        print("\n" + "/home/lyd/webapps/JBrowse/"+srr_n + "\tthe dir exists.\n")
    return [input_dir, ref_dir, output_dir, switch_1, switch_2, switch_3]


def test(srr_n, input_path):
    os.system("mv " + input_path + "/" + srr_n + ".fastq " + input_path + "/" + srr_n + "_download.fastq")
    os.system(
        "tools/seqtk sample -s11 " + input_path + "/" + srr_n + "_download.fastq 10000 > " + input_path + "/" + srr_n + ".fastq")
    os.system("rm " + input_path + "/" + srr_n + "_download.fastq")


def main(argv):
    srr_n = ""
    method = 0
    process_n = 4
    output_path = ""
    kegg_id = ""
    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:m:p:k:",
                                   ["help", "input=", "output=", "method=", "processor=", "keggID="])
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt in ("-i", "--input"):
                srr_n = arg
            elif opt in ("-o", "--output"):
                output_path = arg
            elif opt in ("-m", "--method"):
                method = int(arg)
            elif opt in ("-p", "--processor"):
                process_n = int(arg)
            elif opt in ("-k", "--keggID"):
                kegg_id = arg
        output_path = check_args(srr_n, kegg_id, process_n, method, output_path)
        information = {'Organism': "", 'Instrument': "", 'Layout': "", "keggID": kegg_id}
        getinformation(srr_n, information)
        x = paired_or_single(information["Layout"])
        # print x
        # x = 1 means paired-end and x = 0 means single-end
        _dir = makedir(srr_n, output_path , kegg_id)
        # _dir has 3 dir, input_dir, ref_dir, out_dir
        if _dir[3] == 0:
            sra2fastq(srr_n, x, _dir[0])
            #test(srr_n, _dir[0])
        if _dir[4] == 0:
            download_annotion(kegg_id, _dir[1])
        if _dir[5] == 0:
            if method == 0:
                rockhopper_operon_predict(srr_n, x, _dir, str(process_n))
                # print(_dir[5])
                res2jbrowse(_dir,srr_n,kegg_id)
            elif method == 1:
                # WRITE a python/R/Perl script to connect DOOR and download the .opr file
                download_opr(kegg_id, _dir[1])
                CONDOP_operon_predict(srr_n, x, _dir, str(process_n))
        else:
            print("the result has done !")

    except getopt.GetoptError:
        usage()
        sys.exit(0)


if __name__ == "__main__":
    main(sys.argv)
