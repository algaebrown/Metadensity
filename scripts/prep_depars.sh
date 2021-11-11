# get name
# both ucsc files downloaded from ucsc
ucsc_file=data/hg38_gencode36_wholegene.txt
ucsc_bed=data/hg38_gencode36_wholegene.bed
name_file=${ucsc_file%.txt}.name
cut -d $'\t' -f 4,18  $ucsc_file > $name_file

python src/DaPars_Extract_Anno.py -b $ucsc_bed -s $name_file -o ${ucsc_file}.utr3.bed
