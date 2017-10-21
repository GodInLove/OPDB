import re


def findall_pat(pat, content):
    """
    :method: a method to perfect re.findall()
    :param pat: string
    :param content: string
    :return: list
    """
    pat = re.compile(pat)
    find = re.findall(pat, content)
    return find


def paired_or_single(layout):
    """
    :method: a method to convert paird to 1 and singled to 0
    :param layout: string
    :return: int
    """
    pat_paired = re.compile(r'[padPAD]')
    searchable = re.search(pat_paired, layout)
    if searchable:
        return 1
    else:
        return 0
