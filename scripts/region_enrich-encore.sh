#!/bin/bash
#PBS -q home-yeo
#PBS -N region_enrich_benchmark
#PBS -l nodes=1:ppn=5
#PBS -l walltime=8:00:00
#PBS -o /home/hsher/cluster_msg/region.out
#PBS -e /home/hsher/cluster_msg/region.err
#PBS -t 0-16

ip_bams=(ENCFF391IBH.bam ENCFF501PCN.bam ENCFF923ZPC.bam ENCFF547AHN.bam ENCFF701ZSE.bam ENCFF857OFG.bam ENCFF394THC.bam ENCFF975NJX.bam ENCFF169EKN.bam ENCFF993XXG.bam ENCFF149SFF.bam ENCFF607WID.bam ENCFF341YQM.bam ENCFF583QFB.bam ENCFF177YJA.bam ENCFF701PJM.bam)
input_bams=(ENCFF802LNA.bam ENCFF053BXG.bam ENCFF394DDJ.bam ENCFF190ITO.bam ENCFF062LCG.bam ENCFF949MMS.bam ENCFF405WBU.bam ENCFF810PWF.bam ENCFF615UKX.bam ENCFF209EIX.bam ENCFF878TBN.bam ENCFF255EGS.bam ENCFF447LJU.bam ENCFF222HEX.bam ENCFF117HXE.bam ENCFF696LDT.bam)

ip_bam=${ip_bams[$PBS_ARRAYID]}
input_bam=${input_bams[$PBS_ARRAYID]}

feature=/home/hsher/gencode_coords/gencode.v33.combine.sorted.gff3
bamdir=/home/hsher/seqdata/eclip_raw/
outdir=/home/hsher/seqdata/20200804_charlene_encore_region/

echo $bamdir$ip_bam
echo $bamdir$input_bam
python /home/hsher/projects/Metadensity/metadensity/region_call.py --ip $bamdir$ip_bam --input $bamdir$input_bam --feature $feature --out $outdir${ip_bam%.*}.bed --pval 0.05 --fold 3 --pickle