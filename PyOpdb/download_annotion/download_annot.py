import os

# WRITE a python/R/Perl script to connect NCBI API and download the ref.fna ref.gff ref.rnt ref.gtf
# _dir is a list with (input_path, ref_path, output_path)
import urllib.request

import re

import sys


def callbackfunc(blocknum, blocksize, totalsize):
    percent = 100.0 * blocknum * blocksize / totalsize
    if percent > 100:
        percent = 100
    print("%.2f%%" % percent)


def download_annotion(kegg_id, ref_path):
    print("\ndownloading annotion:.gtf .rnt .gff.......\n")
    # assume that finished it
    # download the annotion in the ref_dir
    url = "http://www.genome.jp/kegg-bin/show_organism?org=" + kegg_id
    print("search in the kegg database\n")
    req = urllib.request.Request(url)
    response = urllib.request.urlopen(req)
    html = response.read()
    html = html.decode('utf-8')
    pat1 = re.compile("(ftp\:\/\/ftp\.ncbi\.nlm\.nih\.gov\/[a-z\/A-Z0-9\_\.]+)")
    result = re.findall(pat1, html)
    if len(result) == 0:
        print("not find the data.please check your input keggID")
    # found the download urls
    name = result[0].split("/")[-1]
    url_download = result[0] + "/" + name
    # fna files download
    print("downloading fna....")
    url_fna = url_download + "_genomic.fna.gz"
    print("\nopening:" + url_fna + "\n")
    local_fna = ref_path + "/" + kegg_id + ".fna.gz"
    urllib.request.urlretrieve(url_fna, local_fna,callbackfunc)
    os.system("gzip -d " + local_fna)
    # gff files download
    print("downloading gff....")
    url_gff = url_download + "_genomic.gff.gz"
    print("\nopening:" + url_gff + "\n")
    local_gff = ref_path + "/" + kegg_id + ".gff.gz"
    urllib.request.urlretrieve(url_gff, local_gff,callbackfunc)
    os.system("gzip -d " + local_gff)
    os.system("mv " + ref_path + "/" + kegg_id + ".gff " + ref_path + "/" + kegg_id + "_orgin.gff")
    # feature files download
    print("downloading feature_table....")
    url_feature = url_download + "_feature_table.txt.gz"
    print("\nopening:" + url_feature + "\n")
    local_feature = ref_path + "/" + kegg_id + "_feature.txt.gz"
    urllib.request.urlretrieve(url_feature, local_feature,callbackfunc)
    os.system("gzip -d " + local_feature)
    os.system("Rscript scripts/prokaryo_anno_download.R " + kegg_id + " " + ref_path)
    print("keggID:" + kegg_id + " , its annotion files were downloaded in the " + ref_path + "\n")
