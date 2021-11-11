#!/bin/bash
#PBS -q hotel
#PBS -N JON_METADEN
#PBS -l nodes=1:ppn=2
#PBS -l walltime=2:50:00
#PBS -o /home/hsher/cluster_msg/jon.out
#PBS -e /home/hsher/cluster_msg/jon.err
#PBS -t 0-0

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metadensity

python ~/Metadensity/scripts/run_metadensity_jon.py /oasis/tscc/scratch/jschmok/2021-09-02_TRNAU1AP_HEK293XT_eCLIP/TRNAU1AP_SE_GRCh38_BothReps/results /oasis/tscc/scratch/jschmok/2021-09-02_TRNAU1AP_HEK293XT_eCLIP/TRNAU1AP_SE_GRCh38_BothReps_MergePeaks/results/01v02.idr.out.normed.bed ~/jonfigures/
