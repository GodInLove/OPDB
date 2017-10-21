import os

import re
import urllib.request

import NCBI_API


def download_annotion(kegg_id, anno_path):
    """
    :method: a method to download annotion files like .fna .gff and .ptt .rnt
    :param kegg_id: string like "eco"
    :param anno_path:
    :return: void(some files in anno_path)
    """
    print("\ndownloading annotion:.gtf .rnt .gff.......\n")
    anno_path = anno_path + kegg_id
    # check whether the anno_path exists or not
    if not os.path.exists(anno_path):
        os.system("mkdir " + anno_path)
    # get the ncbi ftp url through the kegg organisms website
    kegg_url = "http://www.genome.jp/kegg-bin/show_organism?org="
    print("search in the kegg database\n")
    url = kegg_url + kegg_id
    html = NCBI_API.connect_url(url, decode='utf-8')
    url_download = get_url_download(html)
    # downloading fna and gff and convert gff to ptt/rnt
    print("downloading fna....")
    if not os.path.exists(anno_path + "/" + kegg_id + ".fna"):
        # download url
        url_fna = url_download + "_genomic.fna.gz"
        # download path
        local_fna = anno_path + "/" + kegg_id + ".fna.gz"
        # download method
        urllib.request.urlretrieve(url_fna, local_fna, NCBI_API.callback)
    print("downloading gff....")
    if not os.path.exists(anno_path + "/" + kegg_id + "_orgin.gff"):
        url_gff = url_download + "_genomic.gff.gz"
        local_gff = anno_path + "/" + kegg_id + ".gff.gz"
        urllib.request.urlretrieve(url_gff, local_gff, NCBI_API.callback)
        os.system("gzip -d " + local_gff)
        os.system("mv " + anno_path + "/" + kegg_id + ".gff " + anno_path + "/" + kegg_id + "_orgin.gff")
    print("downloading feature_table....")
    if not os.path.exists(anno_path + "/" + kegg_id + ".ptt"):
        url_feature = url_download + "_feature_table.txt.gz"
        local_feature = anno_path + "/" + kegg_id + "_feature.txt.gz"
        urllib.request.urlretrieve(url_feature, local_feature, NCBI_API.callback)
        os.system("gzip -d " + local_feature)
    os.system("Rscript scripts/prokaryo_anno_download.R " + kegg_id + " " + anno_path)
    print("keggID:" + kegg_id + " , its annotion files were downloaded in the " + anno_path + "\n")


def get_url_download(content):
    """
    :method: a sub-method to get the refSeq or GenBank ftp url of a specie through kegg id
    :param content: string
    :return: string
    """
    pat = re.compile("(ftp\:\/\/ftp\.ncbi\.nlm\.nih\.gov\/[a-z\/A-Z0-9\_\.]+)")
    result = re.findall(pat, content)
    name = result[0].split("/")[-1]
    url_download = result[0] + "/" + name
    return url_download
