

genome=/home/hsher/gencode_coords/GRCh38.primary_assembly.genome.fa
/home/hsher/bin/CLAN_release/bin/clan_index -f $genome -d ${genome%.fa}.clan.ref

# bam=$1
# samtools bam2fq -f 4 $bam | sed -n '1~4s/^@/>/p;2~4p' >${bam%.bam}.unmapped.fn
