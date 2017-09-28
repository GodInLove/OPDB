import getopt
import os

import sys

from NCBI_tool.fastq_dump import sra2fastq
from NCBI_tool.getinformation import getinformation, paired_or_single
from check_input.check_argv import usage, check_args
from download_annotion.download_annot import download_annotion
from download_annotion.download_opr import download_opr
from operon_prediction_tool.CONDOP import CONDOP_operon_predict
from operon_prediction_tool.rockhopper import rockhopper_operon_predict


def makedir(srr_n, output_path):
    input_dir = output_path + srr_n + "_input"
    ref_dir = output_path + srr_n + "_ref"
    output_dir = output_path + srr_n + "_output"
    os.system("mkdir " + output_path + srr_n + "_input")
    os.system("mkdir " + output_path + srr_n + "_ref")
    os.system("mkdir " + output_path + srr_n + "_output")
    return [input_dir, ref_dir, output_dir]


def main(argv):
    srr_n = ""
    method = 0
    process_n = 4
    output_path = ""
    nc_n = ""
    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:m:p:n:",
                                   ["help", "input=", "output=", "method=", "processor=", "NC_number="])
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
            elif opt in ("-n", "--NC_number"):
                nc_n = arg
        output_path = check_args(srr_n, nc_n, process_n, method, output_path)
        information = {'Organism': "", 'Instrument': "", 'Layout': "", "NC_number": nc_n}
        getinformation(srr_n, information)
        x = paired_or_single(information["Layout"])
        # print x
        # x = 1 means paired-end and x = 0 means single-end
        _dir = makedir(srr_n, output_path)
        # _dir has 3 dir, input_dir, ref_dir, out_dir
        sra2fastq(srr_n, x, _dir[0])
        download_annotion(srr_n, nc_n, _dir[1])
        if method == 0:
            rockhopper_operon_predict(srr_n, x, _dir, str(process_n))
        elif method == 1:
            # WRITE a python/R/Perl script to connect DOOR and download the .opr file
            download_opr(nc_n, _dir[1])
            CONDOP_operon_predict(srr_n, x, _dir, str(process_n))
        elif method == 2:
            pass
            # RNAseg_operon_predict(srr_n,x,_dir,str(process_n))
    except getopt.GetoptError:
        usage()
        sys.exit(0)


if __name__ == "__main__":
    main(sys.argv)
