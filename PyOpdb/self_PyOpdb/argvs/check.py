import re

import sys

import os

legal = "legal"
illegal = "illegal"


def isLegal(pat, content):
    pat = re.compile(pat)
    searchable = re.search(pat, content)
    if not searchable:
        return illegal
    return legal


def check_output_path(output_path):
    if len(output_path) > 0:
        if output_path[-1] != "/":
            output_path = output_path + "/"
    else:
        print("wrong output_path")
        sys.exit(2)
    if not os.path.exists(output_path):
        os.system("mkdir " + output_path)
    return output_path


def check_srr_n(srr_n):
    return isLegal("[SRA]+[0-9]+", srr_n)


def check_kegg_id(kegg_id):
    return isLegal("[A-Za-z]{3}", kegg_id)


def check_process_n(process_n):
    if 20 > int(process_n) > 0:
        return legal
    return illegal


def check_method(method):
    if 2 > int(method) >= 0:
        return legal
    return illegal


def check_legal_status(legal_status):
    if legal_status == illegal:
        print("please check your input data")
        sys.exit(2)


def check_argvs(srr_n, kegg_id, method, process_n, output_path):
    output_path = check_output_path(output_path)
    legal_status = check_method(method)
    check_legal_status(legal_status)
    legal_status = check_srr_n(srr_n)
    check_legal_status(legal_status)
    legal_status = check_kegg_id(kegg_id)
    check_legal_status(legal_status)
    legal_status = check_process_n(process_n)
    check_legal_status(legal_status)
    return srr_n, kegg_id, method, str(process_n), output_path
