#!/bin/bash
#PBS -q home-yeo
#PBS -N region_motif
#PBS -l nodes=1:ppn=5
#PBS -l walltime=8:00:00
#PBS -o /home/hsher/cluster_msg/homer.out
#PBS -e /home/hsher/cluster_msg/homer.err
#PBS -t 0-25

regions=(ENCFF239CML.bed ENCFF571IRN.bed ENCFF284RSI.bed ENCFF481DSP.bed ENCFF125EPG.bed ENCFF330OFU.bed ENCFF896PKL.bed.intron ENCFF896PKL.bed ENCFF935AZO.bed ENCFF197OBL.bed ENCFF162FHQ.bed ENCFF025YVA.bed ENCFF296GDR.bed ENCFF537RYR.bed ENCFF630YNF.bed ENCFF048BIQ.bed ENCFF646GUE.bed ENCFF068KXE.bed ENCFF830RMG.bed ENCFF896PKL.bed.exon ENCFF170YQV.bed ENCFF876SGL.bed ENCFF756VXR.bed ENCFF690NVV.bed ENCFF492GWE.bed ENCFF233HCD.bed ENCFF959SXJ.bed ENCFF587PLY.bed)
region=${regions[$PBS_ARRAYID]}

outdir=/home/hsher/encore_region_motif/
prefix=${region%.*}
module load homer

# run homer
homer2 denovo -i $outdir$prefix.exon_160.reads.fasta -b $outdir$prefix.exon_input.reads.fasta > $outdir$prefix.exon.motif
homer2 denovo -i $outdir$prefix.intron_160.reads.fasta -b $outdir$prefix.intron_input.reads.fasta > $outdir$prefix.intron.motif