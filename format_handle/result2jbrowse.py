import os


def visual(_dir, kegg_id):
    os.system("cp " + _dir[1] + "/" + kegg_id + ".fna " + _dir[5] + "/" + kegg_id + ".fna")
    os.system("cp " + _dir[1] + "/" + kegg_id + "_orgin.gff " + _dir[5] + "/" + kegg_id + ".gff")
    os.system("rm " + _dir[2] + "/rockhopper/genomeBrowserFiles/" + "_diff*")
    os.system("cp " + _dir[2] + "/rockhopper/genomeBrowserFiles/" + "*.wig " + _dir[5] + "/")


def res2jbrowse(_dir):
    os.system("/home/lyd/Desktop/train_new/web/JBrowse/bin/prepare-refseqs.pl --fasta " + _dir + "/*.fna")
    os.system("/home/lyd/Desktop/train_new/web/JBrowse/bin/flatfile-to-json.pl --gff " + _dir + "/*.gff --trackType CanvasFeatures --trackLabel gff")
    os.system("/home/lyd/Desktop/train_new/web/JBrowse/bin/generate-names.pl -v")
    os.system("for item in " + _dir + "/*wig;do /home/lyd/Desktop/train_new/web/JBrowse/bin/wig-to-json.pl --wig $item --tracklabel ${item##*_};done")