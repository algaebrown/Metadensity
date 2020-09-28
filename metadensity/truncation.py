from pybedtools import BedTool
import pysam


import numpy as np
from scipy.stats import expon, binom
from sklearn.cluster import DBSCAN
import pandas as pd
import math
import sys

from sequence import *


# fetch read start
def read_start_sites(bam_fileobj, interval = None, chrom = None, start = None, end = None, strand = None, single_end = False):
    ''' return read 2 5' end in bedtool interval from bamfile object; strand specific'''
    
    if not single_end:
        if interval:
            subset_reads = list(bam_fileobj.fetch(reference=str(interval.chrom), start=interval.start, end=interval.stop))
            if interval.strand == '+':
                # despite region is enriched with reads, it is possible none of them is read2
                sites = [s.reference_start for s in subset_reads if s.is_read2 and not s.is_reverse]
            else:
                sites = [s.reference_end for s in subset_reads if s.is_read2 and s.is_reverse]
            return sites
        elif chrom:
            subset_reads = list(bam_fileobj.fetch(reference=str(chrom), start=start, end=end))
            if strand == '+':
                sites = [s.reference_start for s in subset_reads if s.is_read2 and not s.is_reverse]
            else:
                sites = [s.reference_end for s in subset_reads if s.is_read2 and s.is_reverse]
    if single_end:
        if interval:
            subset_reads = list(bam_fileobj.fetch(reference=str(interval.chrom), start=interval.start, end=interval.stop))
            if interval.strand == '+':
                sites = [s.reference_start for s in subset_reads if not s.is_reverse]
            else:
                sites = [s.reference_end for s in subset_reads if s.is_reverse]
            return sites
        elif chrom:
            subset_reads = list(bam_fileobj.fetch(reference=str(chrom), start=start, end=end))
            if strand == '+':
                sites = [s.reference_start for s in subset_reads if not s.is_reverse]
            else:
                sites = [s.reference_end for s in subset_reads if s.is_reverse]
    return sites

########################### read start clustering #######################################
class Cluster:
    ''' Read start Cluster'''
    def __init__(self,chrom, start, end, strand, ip_site, input_site):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        
        
        self.ip_site = ip_site
        self.input_site = input_site
        self.no_ip = len(ip_site)
        self.no_input = len(input_site)
        
        self.score = None
        self.seq = ''
    def enrich_score(self, total_ip_in_region, total_input_in_region):
        ''' calculate enrichment score by binom'''
    
    
        # probability of success
        p = total_ip_in_region/(total_ip_in_region + total_input_in_region)
    
        
        # binom
        total_cluster = self.no_ip + self.no_input
    
        bn = binom(total_cluster, p)
        prob = bn.cdf(self.no_ip)
        self.score = prob
        
        # pseudocount
        
        bn = binom(total_cluster+2, p)
        prob = bn.cdf(self.no_ip+1)
        self.pseudoscore = prob
        

    
    
    def to_bedstr(self):
        ''' convert to bedtool compatible string format'''
        bedstr = '{} {} {} {} {} {} {} {}\n'.format(self.chrom, self.start, self.end,0,'.',self.strand, '.', '.')
        return bedstr


def find_cluster(interval, bam_fileobj, inputbam_fileobj, combine = False, eps = 2, min_samples = 2):
    ''' find read start cluster using DBSCAN
    interval: BedTool interval
    bam_fileobj: pysam IP bam
    inputbam_fileobj: pysam Input bam
    combine: if True, create cluster using both IP and Input reads; better for compare to Input.
    eps, min_samples: DBSCAN param
    '''
    # fetch read starts
    sites = read_start_sites(bam_fileobj, interval = interval)
    input_sites = read_start_sites(inputbam_fileobj, interval = interval)
    
    total_ip_in_region = len(sites)
    total_input_in_region = len(input_sites)
    
    if total_ip_in_region == 0:
        # sometimes even region is enriched in coverage, it might contain no 5' sites at all
        #print('site length: {}, input length: {}'.format(len(sites), len(input_sites)))
        return [], input_sites, sites

    # convert to numpy array
    if combine:
        X = np.array(sites + input_sites).reshape(1, -1).T
        identity = np.array(['IP']*len(sites) + ['Input']*len(input_sites)) # remember where each data point comes from
    else:
        X = np.array(sites).reshape(1, -1).T
        identity = np.array(['IP']*len(sites))# remember where each data point comes from

    # run DBSCAN
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(X)

    # convert to Cluster Object
    
    # not all belong to a cluster
    lbl=clustering.labels_[clustering.core_sample_indices_] # labels
    pos=X.flatten()[clustering.core_sample_indices_] # positions
    identity = identity[clustering.core_sample_indices_] # IP or Input
    
    d=pd.DataFrame([lbl, pos, identity], index = ['label', 'position', 'identity']).T
    dmax = d.groupby('label')['position'].max().values.flatten().tolist()
    dmin = d.groupby('label')['position'].min().values.flatten().tolist()
    ip_members = [g[1].loc[g[1]['identity'] == 'IP', 'position'].tolist() for g in d.groupby('label')]
    input_members = [g[1].loc[g[1]['identity'] == 'Input', 'position'].tolist() for g in d.groupby('label')]

    clusters = []
    for start, end, ip_site, input_site in zip(dmin, dmax, ip_members, input_members):
        # create cluster object
        c = Cluster(interval.chrom, start, end, interval.strand, ip_site, input_site)
        
        # calculate enrichment score
        c.enrich_score(total_ip_in_region, total_input_in_region)

        clusters.append(c)
    
    return clusters, input_sites, sites

def control_cluster(clusters, interval):
    ''' given found cluster, return control cluster in the same interval of similar size '''

     # find non-cluster region in interval
    cluster_bedstr = ''

    for c in clusters:
        cluster_bedstr+=c.to_bedstr()
    
    clus_bed = BedTool(cluster_bedstr, from_string = True).saveas()

    left_region = BedTool([interval]).subtract(clus_bed).saveas()

    control_cluster = []  
    lindex = 0
    cindex = 0
    while lindex < len(left_region) and cindex < len(clusters):
        c = clusters[cindex]
        l = left_region[lindex]
        if  l.end-l.start < c.end-c.start:
            # if flanking region is too small
            lindex += 1
        else:
            midpoint = (l.start + l.end)/2
            clus_half_length = (c.end-c.start)/2
            start = math.ceil(midpoint - clus_half_length)
            end = math.ceil(midpoint + clus_half_length)
            
            if end < 0:
                end = 0
            control_cluster.append(Cluster(interval.chrom, start, end, interval.strand, [], [])) # empty input and ip sites
            cindex += 1
            lindex += 1
    
    return control_cluster

