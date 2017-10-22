import sys
import os
import NCBI_API
import self_PyOpdb


def test(srr_n, path):
    """
    :method: a method to test small data with 10000 reads
    :param srr_n: ncbi sra number
    :param path: the path of the sra data
    :return: void(some files in path)
    """
    path = path + srr_n
    os.system("mv " + path + "/" + srr_n + ".fastq " + path + "/" + srr_n + "_download.fastq")
    os.system(
        "tools/seqtk sample -s11 " + path + "/" + srr_n + "_download.fastq 10000 > " + path + "/" + srr_n + ".fastq")
    os.system("rm " + path + "/" + srr_n + "_download.fastq")


if __name__ == "__main__":
    # get argvs form the terminal
    srr_n, kegg_id, method, process_n, output_path = self_PyOpdb.get_argvs(sys.argv)
    # check the argvs to avoid low level error
    srr_n, kegg_id, method, process_n, output_path \
        = self_PyOpdb.check_argvs(srr_n, kegg_id, method, process_n, output_path)
    # srr_information is a dictionary including Organism,Instrument,Layout(singled or paired)
    srr_information = NCBI_API.SRA_API.get_information(srr_n)
    # convert sra data to fastq data with ncbi sra-tookit(fastq-dump)
    self_PyOpdb.srr_n_to_fastq(srr_n, srr_information["Layout"], output_path)
    # test, get 10000 reads from the origin data
    #test(srr_n,output_path)
    # download annotion files like .gff .fna and convert .gff to .ptt|.rnt with ProGenome(a R package)
    NCBI_API.Genome_API.download_annotion(kegg_id, output_path)
    # operon_predict method (code with scripts or function)
    self_PyOpdb.operon_predict(srr_n, kegg_id, srr_information["Layout"], process_n, method, output_path)
