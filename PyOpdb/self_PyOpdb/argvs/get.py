import getopt

import sys


def usage():
    print("\nUseage:python main.py -i SRR5486953 -o /out/path -m 1 -k eco")
    print("\nOr python main.py -i SRR_number -o output_path -m <int> -p <int> -n <string>")
    print("\n-i|--input[string]\tSRR_number(default:null,required)")
    print("\n-o|--output[string]\toutput_dir(default:current dir,not required)")
    print("\n-m|--method[int]\tchoose the software, (0,1,2) means (rockhopper,CONDOP)"
          "(default:0,not required)")
    print("\n-p|--processor[int]\tset the number of processor(default:4,not required)")
    print("\n-k|--keggID[string]\tset the keggID of organism,like eco(default:null,required)")


def get_argvs(argv):
    srr_n = ""
    method = 0
    process_n = 4
    output_path = ""
    kegg_id = ""

    if len(argv) == 1:
        usage()
        sys.exit(0)

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
        return srr_n, kegg_id, method, process_n, output_path
    except getopt.GetoptError:
        usage()
        sys.exit(0)