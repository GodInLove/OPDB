#!/bin/bash

#Usage: sh SRA2OP_0_0.sh SRRXXXX

# Basic V0.0

#1.input SRRXXXX and check it
SRRnumber=$1
# if [SRRnumber==/SRR/]
# then
# 	do
# else
# 	not do
# fi

mkdir ${SRRnumber}_input
mkdir ${SRRnumber}_ref

#2.figure out something and prepare reference data 
'''
figuring out different sequencing platform: Nextseq500/Illumina HiSeq 2500/Illumina HiSeq 2000 and more
figuring out the parameter settings of the software
figuring out whether the data is paired-end or not
figuring out .opr file
figuring out .gff .rnt .fna .gtf
'''
#WRITE a python/R/Perl script to connect NCBI API and get those information

#WRITE a python/R/Perl script to connect NCBI API and download the ref.fna ref.gff ref.rnt ref.gtf

#WRITE a python/R/Perl script to connect DOOR and download the .opr file

#3.download SRR and TO .fastq
./fastq-dump $SRRnumber -O ${SRRnumber}_input/
#./fastq-dump $SRRnumber -split-files -O ${SRRnumber}_input/

#4.considering the quality of the FASTQ data

#5.opting a software artifically or automatically

#example with Rockhopper
#already 
mkdir Rockhopper_${SRRnumber}_output

java -Xmx4g -cp Rockhopper.jar Rockhopper -g ${SRRnumber}_ref/ ${SRRnumber}_input/${SRRnumber}.fastq -o Rockhopper_${SRRnumber}_output/ -TIME
#java -Xmx4g -cp Rockhopper.jar Rockhopper -g ${SRRnumber}_ref/ ${SRRnumber}_input/${SRRnumber}_1.fastq%${SRRnumber}_input/${SRRnumber}_2.fastq -o Rockhopper_${SRRnumber}_output/ -TIME

#6.demonstrating the result
cat Rockhopper_${SRRnumber}_output/*operons.txt

