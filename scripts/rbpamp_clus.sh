#!/bin/bash
#PBS -q home-yeo
#PBS -N rbpamp_clus
#PBS -l nodes=1:ppn=6
#PBS -l walltime=8:00:00
#PBS -o /home/hsher/cluster_msg/rbpampclus.out
#PBS -e /home/hsher/cluster_msg/rbpampclus.err
#PBS -t 0-68

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rbpamp


files=(203_HNRNPC_rep1 204_RBFOX2_rep1 215_TIA1_rep1 216_SRSF9_rep1 218_TRA2A_rep1 240_TRA2A_rep1 247_HNRNPA1_rep1 249_HNRNPU_rep1 281_HNRNPU_rep1 283_HNRNPA1_rep1 285_TIA1_rep1 325_LIN28B_rep1 331_RBM22_rep1 339_TAF15_rep1 340_TARDBP_rep1 345_PCBP2_rep1 387_LIN28B_rep1 388_TAF15_rep1 437_SFPQ_rep1 439_KHSRP_rep1 440_EWSR1_rep1 494_RBM22_rep1 530_RPS3_rep1 540_RPS3_rep1 560_HNRNPL_rep1 571_FUBP3_rep1 614_RPS11_rep1 676_RBFOX2_rep1 678_HNRNPL_rep1 689_FUS_rep1 690_PUM1_rep1 694_FUS_rep1 699_HNRNPC_rep1 755_KHSRP_rep1 203_HNRNPC_rep2 204_RBFOX2_rep2 215_TIA1_rep2 216_SRSF9_rep2 218_TRA2A_rep2 240_TRA2A_rep2 247_HNRNPA1_rep2 249_HNRNPU_rep2 281_HNRNPU_rep2 283_HNRNPA1_rep2 285_TIA1_rep2 325_LIN28B_rep2 331_RBM22_rep2 339_TAF15_rep2 340_TARDBP_rep2 345_PCBP2_rep2 387_LIN28B_rep2 388_TAF15_rep2 437_SFPQ_rep2 439_KHSRP_rep2 440_EWSR1_rep2 494_RBM22_rep2 530_RPS3_rep2 540_RPS3_rep2 560_HNRNPL_rep2 571_FUBP3_rep2 614_RPS11_rep2 676_RBFOX2_rep2 678_HNRNPL_rep2 689_FUS_rep2 690_PUM1_rep2 694_FUS_rep2 699_HNRNPC_rep2 755_KHSRP_rep2)

file=${files[$PBS_ARRAYID]}

python /home/hsher/projects/Metadensity/scripts/rbpamp_clus.py /home/hsher/encore_region_shrink/$file.pickle