import getopt
import os

import sys

from NCBI_tool.fastq_dump import sra2fastq
from NCBI_tool.getinformation import getinformation, paired_or_single
from check_input.check_argv import usage, check_args
from download_annotion.download_annot import download_annotion
from download_annotion.download_opr import download_opr
from operon_prediction_tool.CONDOP import CONDOP_operon_predict
from operon_prediction_tool.RNAseg import RNAseg_operon_predict
from operon_prediction_tool.rockhopper import rockhopper_operon_predict


def makedir(srr_n, output_path):
    status = 0
    input_dir = output_path + srr_n + "_input"
    ref_dir = output_path + srr_n + "_ref"
    output_dir = output_path + srr_n + "_output"
    if not os.path.exists(input_dir):
        os.system("mkdir " + input_dir)
    else:
        status = 1
        print("\nthe file exists.\n")
    if not os.path.exists(ref_dir):
        os.system("mkdir " + ref_dir)
    else:
        status = 1
        print("\nthe file exists.\n")
    if not os.path.exists(output_dir):
        os.system("mkdir " + output_dir)
    else:
        status = 1
        print("\nthe file exists.\n")
    return [input_dir, ref_dir, output_dir, status]


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
        _dir = makedir(srr_n, output_path)
        # _dir has 3 dir, input_dir, ref_dir, out_dir
        if _dir[3] == 0:
            sra2fastq(srr_n, x, _dir[0])
            download_annotion(kegg_id, _dir[1])
        if method == 0:
            os.system("rm "+_dir[1]+"/*.gff")
            rockhopper_operon_predict(srr_n, x, _dir, str(process_n))
        elif method == 1:
            # WRITE a python/R/Perl script to connect DOOR and download the .opr file
            download_opr(kegg_id, _dir[1])
            CONDOP_operon_predict(srr_n, x, _dir, str(process_n))
            pass
        elif method == 2:
            RNAseg_operon_predict(srr_n, x, _dir, str(process_n))
    except getopt.GetoptError:
        usage()
        sys.exit(0)


if __name__ == "__main__":
    main(sys.argv)
