### This whole script is to filter for canonical transcripts. output is gencode_combine_correct

gencode_gff3=$1
prefix=${gencode_gff3%.annotation.gff3}


######################### find exon ####################################
# filter for CDS
awk '{ if ($3 == "exon"){ print } }' $gencode_gff3 > $prefix.exon.gff3
awk '{ if ($3 == "CDS"){ print } }' $gencode_gff3 > $prefix.cds.gff3



######################### find transcript ################################## 

awk '{ if ($3 == "transcript"){ print } }' $gencode_gff3> $prefix.transcript.gff3
# don't remove transcript since it is useful

######################### find intron ################################## 
# transcript - exon = intron (but need to make intron id for yourself
bedtools subtract -a $prefix.transcript.gff3 -b $prefix.exon.gff3 > $prefix.intron.gff3

######################### find UTRs ################################## 
cat $gencode_gff3 | grep prime_UTR > $prefix.utr.gff3
awk '{ if ($3 == "UTR"){ print } }' $gencode_gff3>> $prefix.utr.gff3


######################### COMBINE ################################## 
cat $prefix.cds.gff3 > $prefix.combine.gff3
cat $prefix.exon.gff3 >> $prefix.combine.gff3 # exon includes CDS
cat $prefix.intron.gff3 >> $prefix.combine.gff3
cat $prefix.utr.gff3 >> $prefix.combine.gff3

######################### SORT ################################## 
sort -k1,1 -k4,4n $prefix.combine.gff3 > $prefix.combine.sorted.gff3
# remove intermediate
rm $prefix.combine.gff3


######################### FOR PROTEIN CODING ################################## 
# but these exons/intron include lncRNA and other non-coding RNAs
# use intron_exon_coordinate.sh to filter for protein-coding ones;
# use intron_annotation to get intron ids

# filter for protein_coding only
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $prefix.intron.gff3 > $prefix.proteincode.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $prefix.cds.gff3  >> $prefix.proteincode.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $prefix.utr.gff3  >> $prefix.proteincode.gff3

# sort 
cat $prefix.proteincode.gff3 | sort -k1,1 -k4,4n > $prefix.proteincode.sorted.gff3
rm $prefix.proteincode.gff3

######################## FOR lncRNA ######################################

# filter for lncRNA only
# does not include more complicated biotypes like https://www.gencodegenes.org/pages/biotypes.html
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=lncRNA"){print}}' $prefix.intron.gff3 > $prefix.lncRNA.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=lncRNA"){print}}' $prefix.exon.gff3  >> $prefix.lncRNA.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=lincRNA"){print}}' $prefix.intron.gff3 >> $prefix.lncRNA.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=lincRNA"){print}}' $prefix.exon.gff3  >> $prefix.lncRNA.gff3

# sort 
cat $prefix.lncRNA.gff3 | sort -k1,1 -k4,4n > $prefix.lncRNA.sorted.gff3
rm $prefix.lncRNA.gff3

######################## FOR miRNA ######################################
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=miRNA"){print}}' $prefix.exon.gff3 > $prefix.miRNA.gff3

########################## FOR snRNA ######################################
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=snRNA"){print}}' $prefix.exon.gff3 > $prefix.snRNA.gff3

########################### FOR snoRNA #######################################
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=snoRNA"){print}}' $prefix.exon.gff3 > $prefix.snoRNA.gff3
