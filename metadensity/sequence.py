from pybedtools import BedTool
import pysam
from Bio.Seq import Seq
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

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
        #print('extracting from minus strand')
        
        seq = Seq(seq).reverse_complement().transcribe() # 5' to 3'
    
    else:
        seq = Seq(seq).transcribe()
    
    return seq

############################ kmer related function ##########################################

def kmer(s, k=5):
    ''' given string, return all UNIQUE 5-mer'''
    return list(set([str(s[i:i+k]) for i in range(len(s)) if i+k < len(s)+1]))

def count_kmer(all_seqs, k=5):
    ''' Count k-mer '''
    all_kmer = []
    for s in all_seqs:
        all_kmer.extend(kmer(s, k= k))
    kmer_count = Counter(all_kmer)
    return kmer_count

def kmer_rvalue(ip_seq, input_seq, k= 5):
    ''' Calculate Enrichment Value for each k-mer
    The frequency of kmer in IP / frequency of kmer in Input
    
    '''
    ip_count = count_kmer(ip_seq, k= k)
    input_count = count_kmer(input_seq, k= k)
    
    # number
    #ip_total = sum(ip_count.values())
    #input_total = sum(input_count.values())
    
      
    possible_mer = set(ip_count.keys()).union(input_count.keys())
    
    stat = pd.DataFrame(index = list(possible_mer))
    stat['ip_count'] = stat.index.map(ip_count)
    stat['input_count'] = stat.index.map(input_count)
    stat.fillna(0, inplace = True)
    
    # Calculate frequency
    stat['ip_frequency'] = stat['ip_count']/len(ip_seq)
    stat['input_frequency'] = stat['input_count']/len(input_seq)
    
    stat['R value'] = stat['ip_frequency'].div(stat['input_frequency'])
    
    
    return stat
def enrich_plot(stat, motif, r_thres = 2, freq_thres = 0.125, ax = None):
    ''' plot R value vs IP count kmer enrichment plot '''
    if not ax:
        f, ax = plt.subplots(1,1)
    stat.plot(kind = 'scatter', x = 'R value', y = 'ip_frequency', color = 'grey', ax = ax)
    
    # find high enough points
    high = stat.loc[(stat['R value'] > r_thres) | (stat['ip_frequency'] > freq_thres)]
    for index, row in high.iterrows():
        ax.text(row['R value'], row['ip_frequency'], index)
    
    # highlight canonical motif
    if len(motif) > 0:
        mo = stat.loc[motif]
        mo.plot(kind = 'scatter', x = 'R value', y = 'ip_frequency', color = 'tomato', ax = ax)
        for index, row in mo.iterrows():
            ax.text(row['R value'], row['ip_frequency'], index, color = 'tomato')
    
