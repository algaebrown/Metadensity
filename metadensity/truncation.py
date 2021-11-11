from pybedtools import BedTool
import pysam

from collections import Counter
import numpy as np
from scipy.stats import expon, binom
from sklearn.cluster import DBSCAN
import pandas as pd
import math
import sys

from .sequence import *

def fetch_reads(bam_fileobj, interval = None, chrom = None, start = None, end = None, strand = None, single_end = False, read2 = True):
    ''' fetch reads by pysam.fetch, but strand specific '''
    # fetching reads based on the way we are asking for it
    if interval:
        subset_reads = list(bam_fileobj.fetch(reference=str(interval.chrom), start=interval.start, end=interval.stop))
    elif chrom:
        subset_reads = list(bam_fileobj.fetch(reference=str(chrom), start=start, end=end))
    # paired end, we only fetch read 2
    if single_end:
        pass
    elif read2:
        subset_reads = [s for s in subset_reads if s.is_read2] ### for eCLIP library
    else:
        subset_reads = [s for s in subset_reads if s.is_read1] ### for tao's m6A construction
    # filter for strand
    if strand == '+':
        reads = [s for s in subset_reads if not s.is_reverse]
    else:
        reads = [s for s in subset_reads if s.is_reverse]
    return reads

# fetch read start
def read_start_sites(bam_fileobj, interval = None, chrom = None, start = None, end = None, strand = None, single_end = False, read2 = True):
    ''' return crosslinking events in bedtool interval from bamfile object; strand specific'''
    if single_end:
        # if there is only 1 read, no such read1/read2 issue
        read2=False
    profile = strand_specific_pileup(bam_fileobj, interval = interval, chrom = chrom, start = start, end = end, strand = strand, single_end = single_end, read2 = read2)
    event_count = profile[['trun', 'mismatch', 'del']].sum(axis = 1)
    pos_list = []
    for pos, count in zip(event_count.index.values, event_count.values):
        
        pos_list += [pos]*int(count)
    return pos_list

def truncation_relative_axis(bam_fileobj, interval = None, chrom = None, start = None, end = None, strand = None, single_end = False, read2 = True):
    '''
    return truncation + mismatch + indel count for each position at each site (5' to 3') in np.array()
    '''
    if single_end:
        # if there is only 1 read, no such read1/read2 issue
        read2=False
    if chrom:
        pass
        
    else:
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        strand = interval.strand
    
    # get profile, index = genomic position, columns = number of mismatch, indel and truncation
    profile = strand_specific_pileup(bam_fileobj, chrom=chrom, start=start, end=end, strand=strand, single_end = single_end, read2 = read2)

    site_count = Counter(profile[['trun', 'mismatch', 'del']].sum(axis = 1).to_dict()) # pointing pos -> value #TODO allow selection
        
    pos_count = []
    if strand == '+':
        for pos in range(start, end):
            pos_count.append(site_count[pos])
    if strand == '-':
        for pos in range(end, start, -1): # 15, 14, 13, 12, 11, 10
            pos_count.append(site_count[pos])
    return np.array(pos_count)  
def strand_specific_pileup(bam_fileobj, interval = None, chrom = None, start = None, end = None, strand = None, single_end = False, read2 = True):
    ''' pileup reads of a defined region, 
    return profile (pd.DataFrame), count of truncation, deletion and nucleotides, index = genomic positions '''

    if single_end:
        # if there is only 1 read, no such read1/read2 issue
        read2=False
    
    
    if interval:
        start = interval.start
        end = interval.end
        strand = interval.strand
        chrom = interval.chrom
    
    all_pos_count = []
    all_pos = []
    for pileupcolumn in bam_fileobj.pileup(chrom, start = start, end = end, truncate = True):
    
        pos_profile = []
        for pileupread in pileupcolumn.pileups:
            if (pileupread.alignment.is_reverse and strand == '-' and pileupread.alignment.is_read2 == read2) or ((not pileupread.alignment.is_reverse) and strand == '+' and pileupread.alignment.is_read2 == read2):
                if not pileupread.is_del and not pileupread.is_refskip:
                    # query position is None if is_del or is_refskip is set.
                    pos_profile.append(pileupread.alignment.query_sequence[pileupread.query_position])
            
                if pileupread.indel < 0:
                    pos_profile.append('del')
                    #print(pileupread.indel)
                if pileupread.indel > 0:
                    pos_profile.append('ins')
                if (pileupread.is_tail and strand == '-') or (pileupread.is_head and strand == '+'): # truncation
                    pos_profile.append('trun')
        
        all_pos_count.append(Counter(pos_profile))
        all_pos.append(pileupcolumn.reference_pos)
    profile = pd.DataFrame(all_pos_count, columns = list('ATCG')+['trun', 'del', 'ins'], index = all_pos).fillna(0)
    
    # get mismatch
    get_mismatch(profile)
    return profile

def get_mismatch(profile):
    ''' return # non-majority nucleotide in read '''
    profile['mismatch'] = profile[list('ATCG')].sum(axis = 1)-profile[list('ATCG')].max(axis = 1)
    return profile[list('ATCG')].sum(axis = 1)-profile[list('ATCG')].max(axis = 1)


     


