import re

import sys

import os

legal = "legal"
illegal = "illegal"


def check_argvs(srr_n, kegg_id, method, process_n, output_path):
    """
    :method: a method to check the argvs
    :param srr_n: string
    :param kegg_id: string
    :param method: int
    :param process_n: int
    :param output_path: string
    :return: list
    """
    output_path = check_output_path(output_path)
    print_legal_status(check_method(method), method)
    print_legal_status(check_srr_n(srr_n), srr_n)
    print_legal_status(check_kegg_id(kegg_id), kegg_id)
    print_legal_status(check_process_n(process_n), process_n)
    return srr_n, kegg_id, method, str(process_n), output_path


def isLegal(pat, content):
    """
    :method: a sub-method
    :param pat: string
    :param content: string
    :return: string
    """
    pat = re.compile(pat)
    searchable = re.search(pat, content)
    if not searchable:
        return illegal
    return legal


def check_output_path(output_path):
    """
    :method: a sub-method to make sure that the outpath with a "/"
    :param output_path: string
    :return:  string
    """
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
    """
    :method: a sub-method
    :param srr_n: string
    :return: string
    """
    return isLegal("[SRA]+[0-9]+", srr_n)


def check_kegg_id(kegg_id):
    return isLegal("[A-Za-z]{3,10}", kegg_id)


def check_process_n(process_n):
    if 20 > int(process_n) > 0:
        return legal
    return illegal


def check_method(method):
    if 2 > int(method) >= 0:
        return legal
    return illegal


def print_legal_status(legal_status, data):
    if legal_status == illegal:
        print("please check your input data:" + data)
        sys.exit(2)
