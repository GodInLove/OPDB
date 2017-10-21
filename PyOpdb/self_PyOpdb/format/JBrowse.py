import os


def res2jbrowse(ref_path, result_path, srr_n, kegg_id):
    path = "/home/lyd/webapps/JBrowse"  # the path of Jbrowse
    # handle reference sequence .fasta
    os.system("cd " + path + " && bin/prepare-refseqs.pl --fasta " + ref_path + "/*.fna")
    os.system(
        "cd " + path + " && for item in " + result_path + "/*wig;do bin/wig-to-json.pl --wig $item --tracklabel ${item##*_};done")
    os.system("cp " + result_path + "/*bw " + path + "/data/")
    os.system("cp " + result_path + "/*txt " + path + "/data/")
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_forw_sort.bw --key forv.bw --label rnaseq1 --pos_color blue --plot")
    os.system(
        "cd " + path + " && bin/flatfile-to-json.pl --gff " + ref_path + "/*.gff --trackType CanvasFeatures --trackLabel gff")
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_rev_sort.bw --key rev.bw --label rnaseq2 --pos_color red --plot")
    # os.system(
    #     "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_sort.bw --key all.bw --label rnaseq --pos_color green --neg_color red --plot")
    os.system("cd " + path + " && bin/generate-names.pl -v")
    os.system("cd " + path + " && mv data " + srr_n)
