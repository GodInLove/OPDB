import os


def res2jbrowse(_dir,srr_n):
    os.system("cd /home/lyd/webapps/JBrowse/ && bin/prepare-refseqs.pl --fasta " + _dir[1] + "/*.fna")
    os.system(
        "cd /home/lyd/webapps/JBrowse/ && for item in " + _dir[2] + "/*wig;do bin/wig-to-json.pl --wig $item --tracklabel ${item##*_};done")
    os.system(
        "cd /home/lyd/webapps/JBrowse/ && bin/flatfile-to-json.pl --gff " + _dir[2] + "/*.gff --trackType CanvasFeatures --trackLabel gff")
    os.system("cd /home/lyd/webapps/JBrowse/ && bin/generate-names.pl -v")
    os.system("cd /home/lyd/webapps/JBrowse/ && mv data " + srr_n)
