# filter for protein coding intron/exon/cds
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' gencode_combine_sorted.gff3 > gecode_cds_all.gff3

#awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' gencode_intron.gff3 > gencode_cds_all.gff3

# remain first eight columns
cut -d $'\t' -f 1-8 gecode_cds_all.gff3 > gencode_cds_all_8.gff3

# extract parent transcript id
cut -d $'\t' -f 9 gencode_cds_all.gff3| cut -d ";" -f 4| sed 's/transcript_id=//' > gencode_cds_enst

paste -d "\t" gencode_cds_all_8.gff3 gencode_cds_enst > gencode_cds_all.gff3

rm gencode_cds_all_8.gff3
rm gencode_cds_enst
