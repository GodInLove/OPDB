import os


# WRITE a python/R/Perl script to connect NCBI API and download the ref.fna ref.gff ref.rnt ref.gtf
# _dir is a list with (input_path, ref_path, output_path)
def download_annotion(srr_n, nc_n, ref_path):
    print("downloading annotion:.gtf .rnt .gff.......")
    # assume that finished it
    # download the annotion in the ref_dir
    os.system("cp /home/lyd/Desktop/NC_000913_ref/* " + ref_path)
    print(nc_n + "annotion files were downloaded in the " + ref_path)
