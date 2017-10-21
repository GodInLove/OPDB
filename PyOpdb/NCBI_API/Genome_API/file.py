import os

import re
import urllib.request

import NCBI_API


def get_url_download(content):
    pat = re.compile("(ftp\:\/\/ftp\.ncbi\.nlm\.nih\.gov\/[a-z\/A-Z0-9\_\.]+)")
    result = re.findall(pat, content)
    name = result[0].split("/")[-1]
    url_download = result[0] + "/" + name
    return url_download


def download_annotion(kegg_id, ref_path):
    print("\ndownloading annotion:.gtf .rnt .gff.......\n")
    # ref_path = ref_path + kegg_id
    # if not os.path.exists(ref_path):
    #     os.system("mkdir " + ref_path)
    # kegg_url = "http://www.genome.jp/kegg-bin/show_organism?org="
    # print("search in the kegg database\n")
    # url = kegg_url + kegg_id
    # html = NCBI_API.connect_url(url, decode='utf-8')
    # url_download = get_url_download(html)
    # print("downloading fna....")
    # if not os.path.exists(ref_path + "/" + kegg_id + ".fna"):
    #     url_fna = url_download + "_genomic.fna.gz"
    #     local_fna = ref_path + "/" + kegg_id + ".fna.gz"
    #     urllib.request.urlretrieve(url_fna, local_fna, NCBI_API.callback)
    # print("downloading gff....")
    # if not os.path.exists(ref_path + "/" + kegg_id + "_orgin.gff"):
    #     url_gff = url_download + "_genomic.gff.gz"
    #     local_gff = ref_path + "/" + kegg_id + ".gff.gz"
    #     urllib.request.urlretrieve(url_gff, local_gff, NCBI_API.callback)
    #     os.system("gzip -d " + local_gff)
    #     os.system("mv " + ref_path + "/" + kegg_id + ".gff " + ref_path + "/" + kegg_id + "_orgin.gff")
    # print("downloading feature_table....")
    # if not os.path.exists(ref_path + "/" + kegg_id + ".ptt"):
    #     url_feature = url_download + "_feature_table.txt.gz"
    #     local_feature = ref_path + "/" + kegg_id + "_feature.txt.gz"
    #     urllib.request.urlretrieve(url_feature, local_feature, NCBI_API.callback)
    #     os.system("gzip -d " + local_feature)
    os.system("Rscript scripts/prokaryo_anno_download.R " + kegg_id + " " + ref_path)
    print("keggID:" + kegg_id + " , its annotion files were downloaded in the " + ref_path + "\n")
