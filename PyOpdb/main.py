import sys
import os
import NCBI_API
import self_PyOpdb


def test(srr_n, input_path):
    input_path = input_path + srr_n
    os.system("mv " + input_path + "/" + srr_n + ".fastq " + input_path + "/" + srr_n + "_download.fastq")
    os.system(
        "tools/seqtk sample -s11 " + input_path + "/" + srr_n + "_download.fastq 10000 > " + input_path + "/" + srr_n + ".fastq")
    os.system("rm " + input_path + "/" + srr_n + "_download.fastq")


if __name__ == "__main__":
    srr_n, kegg_id, method, process_n, output_path = self_PyOpdb.get_argvs(sys.argv)
    srr_n, kegg_id, method, process_n, output_path \
        = self_PyOpdb.check_argvs(srr_n, kegg_id, method, process_n, output_path)
    # srr_information is a dictionary including Organism,Instrument,Layout(singled or paired)
    srr_information = NCBI_API.SRA_API.get_information(srr_n)
    self_PyOpdb.srr_n_to_fastq(srr_n, srr_information["Layout"], output_path)
    test(srr_n,output_path)
    NCBI_API.Genome_API.download_annotion(kegg_id, output_path)
    self_PyOpdb.operon_predict(srr_n, kegg_id, srr_information["Layout"], process_n, method, output_path)
