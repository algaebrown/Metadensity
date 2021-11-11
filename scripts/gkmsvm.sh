#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate gkmsvm
#Rscript ~/Metadensity/scripts/gkmsvm.R 

/home/hsher/bin/lsgkm-0.1.1/src/gkmtrain /home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.fasta /home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.ctrl.fasta /home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.test.out -l 7 -k 5 -d 2
/home/hsher/bin/lsgkm-0.1.1/src/gkmpredict /home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.test.fasta /home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.test.out.model.txt /home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.test.output.txt