import os


def res2jbrowse(ref_path, result_path, srr_n, kegg_id):
    path = "/home/lyd/webapps/JBrowse"  # the path of Jbrowse
    # handle reference sequence .fasta
    os.system("cd " + path + " && bin/prepare-refseqs.pl --fasta " + ref_path + "/*.fna")
    # handle result .bw
    os.system("cp " + result_path + "/*bw " + path + "/data/")
    # forv.bw
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_forw_sort.bw --key forv.bw --label rnaseq1 --pos_color blue --plot")
    # rev.bw
    os.system(
        "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_rev_sort.bw --key rev.bw --label rnaseq2 --pos_color red --plot")
    # operon.bw
    os.system(
         "cd " + path + " && bin/add-bw-track.pl --bw_url " + srr_n + "_operon.bw --key all.bw --label rnaseq --pos_color blue --neg_color red --plot")
    # handle reference .gff
    os.system(
        "cd " + path + " && bin/flatfile-to-json.pl --gff " + ref_path + "/*.gff --trackType CanvasFeatures --trackLabel gff")
    # names
    os.system("cd " + path + " && bin/generate-names.pl -v")
    # move rockhopper result
    os.system("cp " + result_path + "/*txt " + path + "/data/")
    # change name "data" to srr_n
    os.system("cd " + path + " && mv data " + srr_n)
