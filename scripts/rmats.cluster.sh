#!/bin/bash
#PBS -q home-yeo
#PBS -N ENCODE_hg38_rmats
#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -o /home/hsher/cluster_msg/encode_rmats.out
#PBS -e /home/hsher/cluster_msg/encode_rmats.err
#PBS -t 0-9


module load rmats/3.2.5
outdir=/projects/ps-yeolab5/encode/rnaseq/alt_splicing_hg38/
basedir=/projects/ps-yeolab5/encode/rnaseq/shrna_knockdown/

bamones=(ENCFF468TYN ENCFF202RTM ENCFF425BAT ENCFF217CYB ENCFF580VRV ENCFF919UOJ ENCFF357UJY ENCFF673GAB ENCFF191YTI ENCFF650ZEJ)

bamtwos=(ENCFF593VAX ENCFF624WAO ENCFF560KLF ENCFF625CKI ENCFF431ARJ ENCFF087WDY ENCFF039VBF ENCFF605AHY nan ENCFF431LMX)

ctrlones=(ENCFF925CTT ENCFF694JWV ENCFF396YDD ENCFF519WJE ENCFF808WAQ ENCFF289TTH ENCFF232HCC ENCFF396YDD ENCFF074QTM ENCFF185IAD)

ctrltwos=(ENCFF360UZH ENCFF902YHV ENCFF867NHV ENCFF272FNP ENCFF181IAL ENCFF284NQQ ENCFF979BGR ENCFF867NHV ENCFF046OMI ENCFF021NQD)

accessions=(ENCSR016IDR ENCSR094KBY ENCSR545AIK ENCSR155EZL ENCSR777EDL ENCSR153GKS ENCSR681SMT ENCSR560AYQ ENCSR829EFL ENCSR490DYI)

lengths=(100 100 100 100 100 100 100 100 100 100)



bam1=${bamones[$PBS_ARRAYID]}
bam2=${bamtwos[$PBS_ARRAYID]}
ctrl1=${ctrlones[$PBS_ARRAYID]}
ctrl2=${ctrltwos[$PBS_ARRAYID]}
acc=${accessions[$PBS_ARRAYID]}
len=${lengths[$PBS_ARRAYID]}

echo $ctrl2

fullout=$outdir$acc
mkdir $fullout

bamall=$basedir$bam1.bam,$basedir$bam2.bam
ctrlall=$basedir$ctrl1.bam,$basedir$ctrl2.bam

rmats \
-b1 $bamall \
-b2 $ctrlall \
-o $fullout \
-t paired \
-len $len \
-gtf /home/hsher/gencode_coords/gencode.v33.annotation.gtf \
-libType fr-firststrand
