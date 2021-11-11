#!/bin/bash
#PBS -q home-yeo
#PBS -N fastq-dump
#PBS -l nodes=1:ppn=6
#PBS -l walltime=20:50:00
#PBS -o /home/hsher/cluster_msg/fastqdump.out
#PBS -e /home/hsher/cluster_msg/fastqdump.err
#PBS -t 0-19

srrs=(SRR11164864 SRR11164865 SRR11164866 SRR11164867 SRR11164868 SRR11164869 SRR11164870 SRR11164871 SRR11164872 SRR11164873 SRR11164874 SRR11164875 SRR11164876 SRR11164877 SRR11164878 SRR11164879 SRR11164880 SRR11164881 SRR11164882 SRR11164883)


srr=${srrs[$PBS_ARRAYID]}

~/bin/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump --outdir /home/hsher/ncbi/fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/sra/$srr.sra
