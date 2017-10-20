import sys
import re
import NCBI_API


def isExist(content, error_print):
    pat_nothing = re.compile(r'No items found')
    searchable = re.search(pat_nothing, content)
    if searchable:
        print(error_print)
        return "no"
    return "yes"


def restrict_match(match):
    if len(match) == 1:
        result = match[0]
        return result
    else:
        print("wrong")
        sys.exit(2)


def restrict_findall(pat, content):
    match = NCBI_API.findall_pat(pat, content)
    return restrict_match(match)


def get_information(srr_n):
    information = {'Organism': "", 'Instrument': "", 'Layout': ""}
    sra_url = "https://www.ncbi.nlm.nih.gov/sra/"

    url = sra_url + "?term=" + srr_n \
          + '%5BAll+Fields%5D+AND+"biomol+rna"%5BProperties%5D'
    print("\nsearch in the SRA database,please wait...\n")
    html = NCBI_API.connect_url(url, decode="utf-8")
    exist_status = isExist(html,
                           error_print="your srr data is not public or not exist")
    if exist_status == "no":
        sys.exit(2)
    pat_org = 'Organism\: \<span\>\<a href\=.*\"\>([a-zA-Z0-9\.\s\-]+).*expand showed sra-full-data'
    information['Organism'] = restrict_findall(pat_org, html)
    pat_ins = 'Instrument\: \<span\>([A-Za-z0-9\s]+)\<\/span\>\<\/div\>\<div\>Strategy\:'
    information['Instrument'] = restrict_findall(pat_ins, html)
    pat_lay = 'Layout: <span>([PAIREDSINGLEpairedsinglenN]+).*sra-full-data'
    information['Layout'] = restrict_findall(pat_lay, html)
    information['Layout'] = NCBI_API.paired_or_single(information['Layout'])
    return information
