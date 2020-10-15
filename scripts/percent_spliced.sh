bam=$1

module load samtools

n_spliced=$(samtools view $bam | awk '($6 ~ /N/)' | cut -f1 | wc -l)
n_total=$(samtools view $bam | cut -f1 | wc -l)


echo $bam $n_spliced $n_total