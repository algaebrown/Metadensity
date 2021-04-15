#!/bin/bash
#PBS -q hotel
#PBS -N EM
#PBS -l nodes=1:ppn=6
#PBS -l walltime=8:00:00
#PBS -o /home/hsher/cluster_msg/em.out
#PBS -e /home/hsher/cluster_msg/em.err
#PBS -t 0-67

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rbpamp

clusters=(247_HNRNPA1_rep2.pickle 203_HNRNPC_rep2.pickle 387_LIN28B_rep1.pickle 345_PCBP2_rep1.pickle 388_TAF15_rep2.pickle 325_LIN28B_rep1.pickle 530_RPS3_rep1.pickle 614_RPS11_rep2.pickle 440_EWSR1_rep1.pickle 755_KHSRP_rep1.pickle 571_FUBP3_rep1.pickle 437_SFPQ_rep1.pickle 215_TIA1_rep2.pickle 249_HNRNPU_rep2.pickle 439_KHSRP_rep2.pickle 218_TRA2A_rep1.pickle 283_HNRNPA1_rep2.pickle 689_FUS_rep1.pickle 690_PUM1_rep1.pickle 285_TIA1_rep1.pickle 216_SRSF9_rep1.pickle 676_RBFOX2_rep1.pickle 494_RBM22_rep1.pickle 694_FUS_rep1.pickle 540_RPS3_rep1.pickle 281_HNRNPU_rep2.pickle 204_RBFOX2_rep2.pickle 340_TARDBP_rep2.pickle 339_TAF15_rep2.pickle 331_RBM22_rep1.pickle 560_HNRNPL_rep2.pickle 240_TRA2A_rep2.pickle 678_HNRNPL_rep1.pickle 699_HNRNPC_rep2.pickle 699_HNRNPC_rep1.pickle 678_HNRNPL_rep2.pickle 240_TRA2A_rep1.pickle 331_RBM22_rep2.pickle 560_HNRNPL_rep1.pickle 340_TARDBP_rep1.pickle 339_TAF15_rep1.pickle 204_RBFOX2_rep1.pickle 281_HNRNPU_rep1.pickle 540_RPS3_rep2.pickle 494_RBM22_rep2.pickle 694_FUS_rep2.pickle 676_RBFOX2_rep2.pickle 216_SRSF9_rep2.pickle 285_TIA1_rep2.pickle 690_PUM1_rep2.pickle 689_FUS_rep2.pickle 283_HNRNPA1_rep1.pickle 218_TRA2A_rep2.pickle 249_HNRNPU_rep1.pickle 439_KHSRP_rep1.pickle 437_SFPQ_rep2.pickle 215_TIA1_rep1.pickle 755_KHSRP_rep2.pickle 571_FUBP3_rep2.pickle 440_EWSR1_rep2.pickle 530_RPS3_rep2.pickle 614_RPS11_rep1.pickle 345_PCBP2_rep2.pickle 388_TAF15_rep1.pickle 325_LIN28B_rep2.pickle 203_HNRNPC_rep1.pickle 387_LIN28B_rep2.pickle 247_HNRNPA1_rep1.pickle)

controls=(247_HNRNPA1_rep2.ctrl.pickle 203_HNRNPC_rep2.ctrl.pickle 387_LIN28B_rep1.ctrl.pickle 345_PCBP2_rep1.ctrl.pickle 388_TAF15_rep2.ctrl.pickle 325_LIN28B_rep1.ctrl.pickle 530_RPS3_rep1.ctrl.pickle 614_RPS11_rep2.ctrl.pickle 440_EWSR1_rep1.ctrl.pickle 755_KHSRP_rep1.ctrl.pickle 571_FUBP3_rep1.ctrl.pickle 437_SFPQ_rep1.ctrl.pickle 215_TIA1_rep2.ctrl.pickle 249_HNRNPU_rep2.ctrl.pickle 439_KHSRP_rep2.ctrl.pickle 218_TRA2A_rep1.ctrl.pickle 283_HNRNPA1_rep2.ctrl.pickle 689_FUS_rep1.ctrl.pickle 690_PUM1_rep1.ctrl.pickle 285_TIA1_rep1.ctrl.pickle 216_SRSF9_rep1.ctrl.pickle 676_RBFOX2_rep1.ctrl.pickle 494_RBM22_rep1.ctrl.pickle 694_FUS_rep1.ctrl.pickle 540_RPS3_rep1.ctrl.pickle 281_HNRNPU_rep2.ctrl.pickle 204_RBFOX2_rep2.ctrl.pickle 340_TARDBP_rep2.ctrl.pickle 339_TAF15_rep2.ctrl.pickle 331_RBM22_rep1.ctrl.pickle 560_HNRNPL_rep2.ctrl.pickle 240_TRA2A_rep2.ctrl.pickle 678_HNRNPL_rep1.ctrl.pickle 699_HNRNPC_rep2.ctrl.pickle 699_HNRNPC_rep1.ctrl.pickle 678_HNRNPL_rep2.ctrl.pickle 240_TRA2A_rep1.ctrl.pickle 331_RBM22_rep2.ctrl.pickle 560_HNRNPL_rep1.ctrl.pickle 340_TARDBP_rep1.ctrl.pickle 339_TAF15_rep1.ctrl.pickle 204_RBFOX2_rep1.ctrl.pickle 281_HNRNPU_rep1.ctrl.pickle 540_RPS3_rep2.ctrl.pickle 494_RBM22_rep2.ctrl.pickle 694_FUS_rep2.ctrl.pickle 676_RBFOX2_rep2.ctrl.pickle 216_SRSF9_rep2.ctrl.pickle 285_TIA1_rep2.ctrl.pickle 690_PUM1_rep2.ctrl.pickle 689_FUS_rep2.ctrl.pickle 283_HNRNPA1_rep1.ctrl.pickle 218_TRA2A_rep2.ctrl.pickle 249_HNRNPU_rep1.ctrl.pickle 439_KHSRP_rep1.ctrl.pickle 437_SFPQ_rep2.ctrl.pickle 215_TIA1_rep1.ctrl.pickle 755_KHSRP_rep2.ctrl.pickle 571_FUBP3_rep2.ctrl.pickle 440_EWSR1_rep2.ctrl.pickle 530_RPS3_rep2.ctrl.pickle 614_RPS11_rep1.ctrl.pickle 345_PCBP2_rep2.ctrl.pickle 388_TAF15_rep1.ctrl.pickle 325_LIN28B_rep2.ctrl.pickle 203_HNRNPC_rep1.ctrl.pickle 387_LIN28B_rep2.ctrl.pickle 247_HNRNPA1_rep1.ctrl.pickle)


clus=${clusters[$PBS_ARRAYID]}
ctrl=${controls[$PBS_ARRAYID]}
inpath=/home/hsher/encore_region_shrink/

# run on 1 high thres 0.5
outpath=/home/hsher/encore_region_em/05
tmpdir=/home/hsher/encore_region_em/05
python /home/hsher/projects/Metadensity/metadensity/em_refine.py --clus $inpath$clus --ctrl $inpath$ctrl --outdir $outpath --tmpdir $tmpdir --logL 0.5

# run on 1 high thres 0.1
outpath=/home/hsher/encore_region_em/01
tmpdir=/home/hsher/encore_region_em/01
python /home/hsher/projects/Metadensity/metadensity/em_refine.py --clus $inpath$clus --ctrl $inpath$ctrl --outdir $outpath --tmpdir $tmpdir --logL 0.1

# run on 1 high thres 0.25
outpath=/home/hsher/encore_region_em/025
tmpdir=/home/hsher/encore_region_em/025
python /home/hsher/projects/Metadensity/metadensity/em_refine.py --clus $inpath$clus --ctrl $inpath$ctrl --outdir $outpath --tmpdir $tmpdir --logL 0.25