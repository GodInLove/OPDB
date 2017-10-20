import re


def findall_pat(pat, content):
    pat = re.compile(pat)
    find = re.findall(pat, content)
    return find

def paired_or_single(layout):
    pat_paired = re.compile(r'[padPAD]')
    searchable = re.search(pat_paired, layout)
    if searchable:
        return 1
    else:
        return 0