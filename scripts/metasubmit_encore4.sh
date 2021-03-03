#!/bin/bash
#PBS -q hotel
#PBS -N metadensity_encode4
#PBS -l nodes=1:ppn=2
#PBS -l walltime=4:50:00
#PBS -o /home/hsher/cluster_msg/meta4.out
#PBS -e /home/hsher/cluster_msg/meta4.err
#PBS -t 0-61

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metadensity

uids=(4001 4002 4004 4007 4008 4009 4006 4010 4012 4014 4017 4018 4019 4020 4022 4023 4035 4037 4038 4040 4044 4047 4052 4053 4053 4050 4059 4045 4054 4055 4058 4056 4048 4049 4098 4099 4110 4030 4036 4082 4093 4072 4081 4028 4083 4087 4097 4096 4111 4107 4108 4109 4114 4116 4117 4068 4065 4084 4088 4104 4126 4070)
uid=${uids[$PBS_ARRAYID]}

#python /home/hsher/projects/Metadensity/scripts/run_metadensity.py $uid ~/seqdata/metadensity

#python /home/hsher/projects/Metadensity/scripts/run_pos_enrich_encode4.py $uid 

#python /home/hsher/projects/Metadensity/scripts/run_kmer_from_read_encode4.py $uid 

#python ~/projects/Metadensity/scripts/run_shape_from_read_encode4.py $uid

#python ~/Metadensity/scripts/rg4_enrichment.py $uid

python ~/Metadensity/scripts/run_metadensity.py $uid ~/densities/

