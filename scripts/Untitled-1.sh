#!/bin/bash
#PBS -q home-yeo
#PBS -N tao_region_enrich
#PBS -l nodes=1:ppn=5
#PBS -l walltime=100:00:00
#PBS -o /home/hsher/cluster_msg/tao_region.out
#PBS -e /home/hsher/cluster_msg/tao_region.err
#PBS -t 0-2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metadensity

ip_bams = (TGIRTAligned.sortedByCoord.out.bam SSIIIAligned.sortedByCoord.out.bam adaptor2Aligned.sortedByCoord.out.bam)
input_bam=InputAligned.sortedByCoord.out.bam

ip_bam=${ip_bams[$PBS_ARRAYID]}


feature=/home/hsher/gencode_coords/gencode.v33.combine.sorted.gff3
bamdir=/home/hsher/tao_sc_m6A/complete_bam/
outdir=/home/hsher/tao_sc_m6A/complete_bam/region_call/

echo $bamdir$ip_bam
echo $bamdir$input_bam
python /home/hsher/projects/Metadensity/metadensity/region_call.py --ip $bamdir$ip_bam --input $bamdir$input_bam --feature $feature --out $outdir${ip_bam%.*}.bed