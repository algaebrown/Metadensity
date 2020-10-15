peak=$1
bed=${peak##*/}
prefix=${bed%.*}
coords=$2 # prefix like gencode.v33 or gencode.vM25
genome_fasta=$3
outpath=$4

coordpath=~/gencode_coords/

mkdir $outpath
# split peaks into intron, exon, UTRs
grep "exon" $peak > $outpath$prefix.exon
grep "UTR3" $peak > $outpath$prefix.UTR3
grep "UTR5" $peak > $outpath$prefix.UTR5
grep -v "exon" $peak | grep -v "UTR3" | grep -v "UTR5" > $outpath$prefix.intron

# get input
bedtools intersect -a $coordpath$coords.exon.gff3 -b $outpath$prefix.exon -v > $outpath$prefix.exon.input
bedtools intersect -a $coordpath$coords.utr.gff3 -b $outpath$prefix.UTR3 -v > $outpath$prefix.UTR3.input
bedtools intersect -a $coordpath$coords.utr.gff3 -b $outpath$prefix.UTR5 -v > $outpath$prefix.UTR5.input
bedtools intersect -a $coordpath$coords.intron.gff3 -b $outpath$prefix.intron -v > $outpath$prefix.intron.input

# extract sequence

for f in $outpath$prefix*
do
    if echo "$f" | grep "input"
    then
        outf=${f%.*}_input.reads
    else
        outf=${f}_160.reads
    fi
echo $outf
echo "filename"
bedtools getfasta -bed $f -fi $coordpath$genome_fasta -s > $outf.fasta
grep -v ">" $outf.fasta > $outf
done

