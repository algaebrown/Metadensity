module load samtools
bam=$1
outprefix=$2
samtools view -h -s 5.75 $bam | samtools view -Sb - > $outprefix.75.bam
samtools view -h -s 5.50 $bam | samtools view -Sb - > $outprefix.50.bam
samtools view -h -s 5.25 $bam | samtools view -Sb - > $outprefix.25.bam

samtools index $outprefix.75.bam
samtools index $outprefix.25.bam
samtools index $outprefix.50.bam
