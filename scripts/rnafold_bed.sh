bed_or_pickle=$1
outdir=$2

mkdir $outdir
rm $outdir/*

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metadensity

python /home/hsher/Metadensity/scripts/prepare_rnafold.py $bed_or_pickle $outdir

conda deactivate

conda activate rnatools
cd $outdir
for f in *.shape.dat
do
    RNAfold --shape $f < ${f%.shape.dat}.fasta > ${f%.shape.dat}.fold
    python /home/hsher/Metadensity/scripts/plot_csln_on_rna.py ${f%.shape.dat}.fold ${f%.shape.dat}.csv
done

