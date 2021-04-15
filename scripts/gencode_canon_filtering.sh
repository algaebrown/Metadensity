### This whole script is to filter for canonical transcripts. output is gencode_combine_correct

gencode_gff3=$1
prefix=${gencode_gff3%.annotation.gff3}
known_canon=$2

######################### find exon ####################################
# filter for canonical transcript (knownCanonical) in the whole gencode
awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ":"); if(fields[2] in save) {print}}' $known_canon $gencode_gff3 > $prefix.filtered.exon.gff3
awk '{ if ($3 == "exon"){ print } }' $prefix.filtered.exon.gff3 > $prefix.exon.gff3

awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ";"); a=fields[2]; b=gensub(/Parent=/, "", 1, a);if(b in save) {print}}' $known_canon $gencode_gff3 > $prefix.filtered.cds.gff3
awk '{ if ($3 == "CDS"){ print } }' $prefix.filtered.cds.gff3 > $prefix.cds.gff3

rm $prefix.filtered.exon.gff3
rm $prefix.filtered.cds.gff3

######################### find transcript ################################## 
awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ";"); a = fields[1];b = gensub(/ID=/, "", 1, a);if(b in save) {print}}' $known_canon $gencode_gff3 > $prefix.filtered.transcript.gff3
awk '{ if ($3 == "transcript"){ print } }' $prefix.filtered.transcript.gff3> $prefix.transcript.gff3
# don't remove transcript since it is useful

######################### find intron ################################## 
# transcript - exon = intron (but need to make intron id for yourself
bedtools subtract -a $prefix.transcript.gff3 -b $prefix.exon.gff3 > $prefix.intron.gff3

######################### find UTRs ################################## 
cat $gencode_gff3 | grep prime_UTR > $prefix.unfiltered.utr.gff3
awk '{ if ($3 == "UTR"){ print } }' $gencode_gff3>> $prefix.unfiltered.utr.gff3 # for gencode v19
awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ";"); a = fields[4];b = gensub(/transcript_id=/, "", 1, a) ;if(b in save) {print}}' $known_canon $prefix.unfiltered.utr.gff3 > $prefix.utr.gff3



# remove intermediate
rm $prefix.unfiltered.utr.gff3

######################### COMBINE ################################## 
cat $prefix.exon.gff3 > $prefix.combine.gff3
cat $prefix.cds.gff3 >> $prefix.combine.gff3
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
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $prefix.exon.gff3  >> $prefix.proteincode.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $prefix.utr.gff3  >> $prefix.proteincode.gff3

# sort 
cat $prefix.proteincode.gff3 | sort -k1,1 -k4,4n > $prefix.cds.sorted.gff3
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
