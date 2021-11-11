module load samtools
bam=$1
outprefix=$2

# calculate the fraction
mapped=$(samtools view -c -F 260 $bam)
target=2000000

morereads=$(echo "$mapped > $target" | bc -l)
if (($morereads))
    then 
    frac=$(echo "scale=2; {$target/$mapped}" | bc )
    samtools view -h -s 5$frac $bam | samtools view -Sb - > $outprefix.twomillion.bam
    echo "downsample to 2 million"

else
    # directly copy without sampling
    #cp $bam $outprefix.twomillion.bam
    echo "WARNING: less than 2 million reads"

fi


samtools index $outprefix.twomillion.bam

