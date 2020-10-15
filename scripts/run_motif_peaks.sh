#!/bin/bash
#PBS -q home-yeo
#PBS -N region_motif
#PBS -l nodes=1:ppn=5
#PBS -l walltime=8:00:00
#PBS -o /home/hsher/cluster_msg/motif.out
#PBS -e /home/hsher/cluster_msg/motif.err
#PBS -t 0-25


regions=(ENCFF239CML.bed ENCFF571IRN.bed ENCFF284RSI.bed ENCFF481DSP.bed ENCFF125EPG.bed ENCFF330OFU.bed testENCFF896PKL.bed.intron ENCFF896PKL.bed ENCFF935AZO.bed ENCFF197OBL.bed ENCFF162FHQ.bed ENCFF025YVA.bed ENCFF296GDR.bed ENCFF537RYR.bed ENCFF630YNF.bed ENCFF048BIQ.bed ENCFF646GUE.bed ENCFF068KXE.bed ENCFF830RMG.bed testENCFF896PKL.bed.exon ENCFF170YQV.bed ENCFF876SGL.bed ENCFF756VXR.bed ENCFF690NVV.bed ENCFF492GWE.bed ENCFF233HCD.bed ENCFF959SXJ.bed ENCFF587PLY.bed)
region=${regions[$PBS_ARRAYID]}

indir=/home/hsher/encore_region_call/
outdir=/home/hsher/encore_region_motif/

# create fasta for homer
#~/projects/Metadensity/scripts/motif.sh  ENCFF896PKL.bed gencode.v33 GRCh38.p13.genome.fa test/
prefix=${region%.*}
outpath=$outdir$prefix/
/home/hsher/projects/Metadensity/scripts/motif.sh  $indir$region gencode.v33 GRCh38.p13.genome.fa $outpath


# run homer
homer2 denovo -i $outpath$prefix.exon_160.reads.fasta -b $outpath$prefix.exon_input.fasta > $outpath$prefix.exon.motif
homer2 denovo -i $outpath$prefix.intron_160.reads.fasta -b $outpath$prefix.intron_input.fasta > $outpath$prefix.intron.motif