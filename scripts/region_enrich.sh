#!/bin/bash
#PBS -q home-yeo
#PBS -N YTHDF_region_enrich
#PBS -l nodes=1:ppn=5
#PBS -l walltime=100:00:00
#PBS -o /home/hsher/cluster_msg/ythdf_region.out
#PBS -e /home/hsher/cluster_msg/ythdf_region.err
#PBS -t 0-11

ip_bams=(zt2_liver_YTHDF2_eCLIP.m1_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m2_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m4_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m11_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m12_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m9_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m6_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m5_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m3_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m10_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m8_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m7_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam)
input_bams=(zt2_liver_YTHDF2_eCLIP.m1_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m2_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m4_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m11_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m12_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m9_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m6_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m5_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_YTHDF2_eCLIP.m3_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m10_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m8_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam zt2_liver_m6A_eCLIP.m7_Input.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam)

ip_bam=${ip_bams[$PBS_ARRAYID]}
input_bam=${input_bams[$PBS_ARRAYID]}

feature=/home/hsher/gencode_coords/gencode.vM25.combine.sorted.gff3
bamdir=/home/hsher/YTHDF2/bams/
outdir=/home/hsher/YTHDF2/region_call/

echo $bamdir$ip_bam
echo $bamdir$input_bam
python /home/hsher/projects/Metadensity/metadensity/region_call.py --ip $bamdir$ip_bam --input $bamdir$input_bam --feature $feature --out $outdir${ip_bam%.*}.bed