
from .truncation import read_start_sites
from .sequence import get_truncation_seq, simulate_kmer_background, kmer_zscore
import pandas as pd

def get_all_site_seq(bam, chrom = 'chr1', start = 0, end = 248956422, strand = '+', window = 25,  single_end = False, read2 = True):
    ''' fetch sequence around read start sites'''
    seq_around = []
    
    sites = read_start_sites(bam, chrom = chrom, start = start, end = end, strand = strand, single_end = single_end, read2 = read2)
    seq = [get_truncation_seq(chrom, s, strand, window = window) for s in sites]

    return seq

def main(bam_rep1, bam_rep2, bam_input1, bam_input2 = None, k=7, chrom = 'chr1', start = 0, end = 248956422, strand = '+', window = 25, n_sample = 1000, n_iter = 100, single_end = False):
    # IPs
    print('fetching IP sequences')
    rep1_seqs = get_all_site_seq(bam_rep1, chrom = chrom, start = start, end = end, strand = strand, window = window, single_end = single_end)
    rep2_seqs = get_all_site_seq(bam_rep2, chrom = chrom, start = start, end = end, strand = strand, window = window, single_end = single_end)

    # Inputs
    print('fetching Input sequences')
    if bam_input2:
        # 2 IP, 2 Input ENCODE 4 structure
        input1_seqs = get_all_site_seq(bam_input1, chrom = chrom, start = start, end = end, strand = strand, window = window, single_end = single_end)
        input2_seqs = get_all_site_seq(bam_input2, chrom = chrom, start = start, end = end, strand = strand, window = window, single_end = single_end)
    else:
        input1_seqs = get_all_site_seq(bam_input1, chrom = chrom, start = start, end = end, strand = strand, window = window, single_end = single_end)
    print('get {}, {} sequences'.format(len(rep1_seqs), len(rep2_seqs)))

    # run seperatly for 2 reps, and combined, compare correlation
    print('Simulating background k-mer from Input')
    bg1 = simulate_kmer_background(input1_seqs, k= k, n_sample = n_sample, n_iter = n_iter)
    if bam_input2:
        bg2 = simulate_kmer_background(input2_seqs, k= k, n_sample = n_sample, n_iter = n_iter)
        bg_combine = simulate_kmer_background(input1_seqs + input2_seqs, k= k, n_sample = n_sample, n_iter = n_iter)
    
    # get z-score
    print('Calculating Z-score')
    if bam_input2:
        z1=kmer_zscore(rep1_seqs, bg1[0], bg1[1], k = k) # mean and std ing bg
        z2=kmer_zscore(rep2_seqs, bg2[0], bg2[1], k = k)
        z_combine = kmer_zscore(rep1_seqs+rep2_seqs, bg_combine[0], bg_combine[1], k = k)
    else:
        z1=kmer_zscore(rep1_seqs, bg1[0], bg1[1], k = k) # mean and std ing bg
        z2=kmer_zscore(rep2_seqs, bg1[0], bg1[1], k = k)
        z_combine = kmer_zscore(rep1_seqs+rep2_seqs, bg1[0], bg1[1], k = k)
    print(len(z1), len(z2), len(z_combine))

    return pd.concat([z1, z2, z_combine], axis = 1)




