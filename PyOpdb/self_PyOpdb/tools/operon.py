import os

from self_PyOpdb.format import res2jbrowse, extract_Synonym, extract_wig
from self_PyOpdb.tools.sam import samtools, bamCoverage


def rockhopper_operon_predict(srr_n, layout, kegg_id, process_n, output_path):
    _dir = [output_path + srr_n, output_path + kegg_id, output_path + "rockhopper" + srr_n]
    if not os.path.exists(_dir[2]):
        os.system("mkdir " + _dir[2])
        print("running Rockhopper....")
        # if it is paired
        if layout == 1:
            os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
                0] + "/" + srr_n + "_1.fastq%" + _dir[0] + "/" + srr_n + "_2.fastq -o " + _dir[
                          2] + " -TIME -SAM")
        else:
            os.system("java -Xmx3g -cp tools/Rockhopper.jar Rockhopper -p " + process_n + " -g " + _dir[1] + " " + _dir[
                0] + "/" + srr_n + ".fastq -o " + _dir[2] + " -TIME -SAM")
        # handle some output files to unify format
        extract_Synonym(_dir[1],_dir[2])
        os.system("rm " + _dir[2] + "/genomeBrowserFiles/*_diff*")
        os.system("mv " + _dir[2] + "/genomeBrowserFiles/*_operons.wig " + _dir[2])
        gff_path = _dir[1] + "/" + kegg_id + ".gff"
        if os.path.exists(_dir[1] + "/" + kegg_id + "_orgin.gff"):
            os.system("rm " + gff_path)
            os.system("cp " + _dir[1] + "/" + kegg_id + "_orgin.gff " + gff_path)
        # convert wig to bigwig
        extract_wig(_dir[2]+"/_operons.wig", gff_path, _dir[2])
    else:
        print("\nthe result has done ! Please check the path： " + _dir[2] + "\n")


def operon_predict(srr_n, kegg_id, layout, process_n, method, output_path):
    # rockhopper
    if method == 0:
        rockhopper_operon_predict(srr_n, layout, kegg_id, process_n, output_path)
        samtools(srr_n, output_path+"rockhopper"+srr_n)
        bamCoverage(srr_n, output_path+"rockhopper"+srr_n)
        res2jbrowse( output_path + kegg_id, output_path + "rockhopper"+srr_n
                    ,srr_n,kegg_id)
    # CONDOP
    elif method == 1:
        pass
    else:
        print("wrong method choose:"+ method+"\nplease use 0 or 1")



