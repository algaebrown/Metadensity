from pybedtools import BedTool
import pysam
from Bio.Seq import Seq

# human genome here
genome='/home/hsher/gencode_coords/GRCh38.p13.genome.fa'
fasta = pysam.FastaFile(genome)

# basic function to fetch sequence strand specific ways

def getRNAsequence(interval, fasta = fasta):
    ''' fetch RNA sequence from bedtool interval'''
    seq = fasta.fetch(reference = interval.chrom, start=interval.start, end=interval.stop)
    
    if interval.strand == '-':
        
        seq = Seq(seq).reverse_complement().transcribe()
    
    else:
        seq = Seq(seq).transcribe()
    
    return seq

def get_truncation_seq(chrom, start, strand, window = 10 , fasta = fasta):
    ''' get sequence around window of truncation site'''
    seq = fasta.fetch(reference = chrom, start=start - window, end = start+window)
    if strand == '-':
        
        seq = Seq(seq).reverse_complement().transcribe()
    
    else:
        
        seq = Seq(seq).transcribe()
    
    return seq
def get_interval_seq(chrom, start, end, strand, fasta = fasta):
    ''' get sequence around interval'''
    seq = fasta.fetch(reference = chrom, start=start, end=end)
    
    if strand == '-':
        
        seq = Seq(seq).reverse_complement().transcribe() # 5' to 3'
    
    else:
        seq = Seq(seq).transcribe()
    
    return seq