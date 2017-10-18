import os


def res2jbrowse(_dir,srr_n,kegg_id):
    os.system("cd /home/lyd/webapps/JBrowse/ && bin/prepare-refseqs.pl --fasta " + _dir[1] + "/*.fna")
    os.system(
        "cd /home/lyd/webapps/JBrowse/ && for item in " + _dir[2] + "/*wig;do bin/wig-to-json.pl --wig $item --tracklabel ${item##*_};done")
    if os.path.exists(_dir[1] + "/" + kegg_id + "_orgin.gff"):
        os.system("rm " + _dir[1] +"/" + kegg_id +".gff")
        os.system("cp " + _dir[1] + "/" + kegg_id + "_orgin.gff " + _dir[1] +"/" + kegg_id +".gff")
    os.system(
        "cd /home/lyd/webapps/JBrowse/ && bin/flatfile-to-json.pl --gff " + _dir[1] + "/*.gff --trackType CanvasFeatures --trackLabel gff")
    os.system("mv " + _dir[2] + "/*bw /home/lyd/webapps/JBrowse/data/")
    os.system("cd /home/lyd/webapps/JBrowse/ && bin/add-bw-track.pl --bw_url " + srr_n + "_forw_sort.bw --key forv.bw --label rnaseq")
    os.system(
        "cd /home/lyd/webapps/JBrowse/ && bin/add-bw-track.pl --bw_url " + srr_n + "_rev_sort.bw --key rev.bw --label rnaseq")
    os.system("cd /home/lyd/webapps/JBrowse/ && bin/generate-names.pl -v")
    os.system("cd /home/lyd/webapps/JBrowse/ && mv data " + srr_n)

