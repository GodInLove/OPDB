import re

import sys


def usage():
    print("\nUseage:python SRA2OP_0_0.py -i SRR5486953 -o /home/lyd/Desktop -m 1 -k eco")
    print("\nOr python SRA2OP_0_0.py -i SRR_number -o output_path -m <int> -p <int> -n <string>")
    print("\n-i|--input[string]\tSRRnumber(default:null)")
    print("\n-o|--output[string]\toutput_dir(default:current dir)")
    print("\n-m|--method[int]\tchoose the software, (0,1,2) means (rockhopper,CONDOP,RNAseg)(default:0)")
    print("\n-p|--processor[int]\tset the number of processor(default:4)")
    print("\n-k|--keggID[string]\tset the keggID of organism,like eco(default:null)")


def checksrr_n(srr_n):
    pat = re.compile(r'[SRA]+[0-9]+')
    searchable = re.search(pat, srr_n)
    if not searchable:
        print("wrong SRA_number!")
        usage()
        sys.exit(2)


def checknc_n(kegg_id):
    pat = re.compile(r'[A-Za-z]{3}')
    searchobj = re.search(pat, kegg_id)
    if not searchobj:
        print("wrong kegg_id!")
        usage()
        sys.exit(2)


# def check_SRR(information):
# 	#check whether the SRR data is up to standard or not
# 	if information['Organism'] == "Escherichia coli":
# 		not_qualified()
# 	else:
# 		pass

def checkprocess_n(process_n):
    if int(process_n) < 0:
        print("wrong process number!")
        usage()
        sys.exit(2)
    elif int(process_n) > 100:
        print("wrong process number!")
        usage()
        sys.exit(2)


def checkmethod(method):
    if int(method) < 0:
        print("wrong method number!")
        usage()
        sys.exit(2)
    elif int(method) > 2:
        print("wrong method number!")
        usage()
        sys.exit(2)


def checkoutput_path(output_path):
    if len(output_path) > 0:
        if output_path[-1] != "/":
            output_path = output_path + "/"
    return output_path


def check_args(srr_n, kegg_id, process_n, method, output_path):
    checkmethod(method)
    checknc_n(kegg_id)
    checkprocess_n(process_n)
    output_path = checkoutput_path(output_path)
    checksrr_n(srr_n)
    return output_path
