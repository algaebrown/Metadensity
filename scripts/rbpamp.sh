
#!/bin/bash
#PBS -q home-yeo
#PBS -N rbpamp_kmer
#PBS -l nodes=1:ppn=2
#PBS -l walltime=8:00:00
#PBS -o /home/hsher/cluster_msg/rbpamp_kmer.out
#PBS -e /home/hsher/cluster_msg/rbpamp_kmer.err
#PBS -t 0-1


source ~/miniconda3/etc/profile.d/conda.sh
conda activate rbpamp

python /home/hsher/projects/Metadensity/scripts/rbpamp_kmer.py