import re
import urllib.request
import sys


def download_opr(kegg_id, ref_path):
    url_kegg = "http://www.genome.jp/kegg-bin/show_organism?org=" + kegg_id
    print("search in the kegg database\n")
    req = urllib.request.Request(url_kegg)
    response = urllib.request.urlopen(req)
    html = response.read()
    html = html.decode('utf-8')
    pat = re.compile("\>(NC\_[0-9]+)\<")
    result_kegg = re.findall(pat,html)
    if len(result_kegg) == 0:
        print("the keggID has not a NC_number")
        sys.exit(2)
    nc_n = result_kegg[0]
    root_url = "http://csbl.bmb.uga.edu/DOOR/"
    # browser search string
    try:
        url = "http://csbl.bmb.uga.edu/DOOR/search_ajax.php?keyword=" + nc_n + "&mode=DataTable&sEcho=1&iColumns=6&sColumns=&iDisplayStart=0&iDisplayLength=15&mDataProp_0=0&mDataProp_1=1&mDataProp_2=2&mDataProp_3=3&mDataProp_4=4&mDataProp_5=5&_=1506424463060"
        print("search in the DOOR database,please wait...\n")
        # print(url)
        req = urllib.request.Request(url)
        response = urllib.request.urlopen(req)
        html = response.read()
        # print (html)
        # check whether the opr file is available or not
        pat_nothing = re.compile(r'iTotalRecords\"\:([0-9]+)')
        html = html.decode('gbk')
        searchable = re.findall(pat_nothing, html)
        if searchable[0] == '0':
            print("the DOOR database has not relative .opr file,please change another method")
            sys.exit(2)
        # when the opr file is available
        _str = html.split("[")
        pat1 = re.compile(r'(operon\.php\?id\=[0-9]+).*' + nc_n)
        for item in _str:
            result1 = re.findall(pat1, item)
            if result1:
                break
        # print(result1)
        url1 = root_url + result1[0]
        print("search in the DOOR database,please wait...\n")
        req = urllib.request.Request(url1)
        response = urllib.request.urlopen(req)
        html = response.read()
        html = html.decode('gbk')
        pat2 = re.compile(r'\<li\>\<a\shref\=\"(.*)\".*\<\/li\>\s+\<li\>\<a\shref\=\"operon\.php\?id\=')
        result2 = re.findall(pat2, html)
        if len(result2) == 0:
            print("cannot find the url,please contact with the author")
        url2 = root_url + result2[0]
        print("search in the DOOR database,please wait...\n")
        req = urllib.request.Request(url2)
        response = urllib.request.urlopen(req)
        html = response.read()
        html = html.decode('gbk')
        pat3 = re.compile(r'window\.location\=\'(downloadNCoperon.php.*)\'\"\svalue\=\"Download\soperon\stable\"\>')
        result3 = re.findall(pat3, html)
        if len(result3) == 0:
            print("cannot find the url,please contact with the author")
        url_download = root_url + result3[0]
        pat_nc = re.compile(r'[0-9]+')
        result = re.findall(pat_nc, url_download)
        local = ref_path + "/" + result[0] + ".opr"
        print("starting download the .opr file,please wait...\n")
        urllib.request.urlretrieve(url_download, local)
        print("NC_number:" + nc_n + " , its opr files were downloaded in the " + ref_path+ "\n")
    except:
        print("wrong!!",sys.exc_info()[0])
        sys.exit(2)
