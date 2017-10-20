import os


def res2jbrowse(_dir, srr_n, kegg_id):
    path = "/home/lyd/webapps/JBrowse"
    os.system("cd " + path + " && bin/prepare-refseqs.pl --fasta " + _dir[1] + "/*.fna")
    os.system(
        "cd " + path + " && for item in " + _dir[
            2] + "/*wig;do bin/wig-to-json.pl --wig $item --tracklabel ${item##*_};done")
    if os.path.exists(_dir[1] + "/" + kegg_id + "_orgin.gff"):
        os.system("rm " + _dir[1] + "/" + kegg_id + ".gff")
        os.system("cp " + _dir[1] + "/" + kegg_id + "_orgin.gff " + _dir[1] + "/" + kegg_id + ".gff")
    os.system("cp " + _dir[2] + "/*bw " + path + "/data/")
    os.system("cp " + _dir[2] + "/*txt " + path + "/data/")
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_forw_sort.bw --key forv.bw --label rnaseq1 --pos_color blue --plot")
    os.system(
        "cd " + path + " && bin/flatfile-to-json.pl --gff " + _dir[
            1] + "/*.gff --trackType CanvasFeatures --trackLabel gff")
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_rev_sort.bw --key rev.bw --label rnaseq2 --pos_color red --plot")
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_sort.bw --key all.bw --label rnaseq --pos_color green --neg_color red --plot")
    os.system("cd " + path + " && bin/generate-names.pl -v")
    os.system("cd " + path + " && mv data " + srr_n)
