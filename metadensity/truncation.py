from pybedtools import BedTool
import pysam


import numpy as np
from scipy.stats import expon, binom
from sklearn.cluster import DBSCAN
import pandas as pd
import math
import sys
sys.path.append('/home/hsher/projects/Metadensity')
from metadensity.sequence import *


# fetch read start
def read_start_sites(bam_fileobj, interval = None, chrom = None, start = None, end = None, strand = None):
    ''' return read 2 5' end in bedtool interval from bamfile object; strand specific'''
    if interval:
        subset_reads = list(bam_fileobj.fetch(reference=str(interval.chrom), start=interval.start, end=interval.stop))
        if interval.strand == '+':
            sites = [s.reference_start for s in subset_reads if s.is_read2]
        else:
            sites = [s.reference_end for s in subset_reads if s.is_read2]
        return sites
    elif chrom:
        subset_reads = list(bam_fileobj.fetch(reference=str(chrom), start=start, end=end))
        if strand == '+':
            sites = [s.reference_start for s in subset_reads if s.is_read2]
        else:
            sites = [s.reference_end for s in subset_reads if s.is_read2]
        return sites

########################### read start clustering #######################################
class Cluster:
    ''' Read start Cluster'''
    def __init__(self,chrom, start, end, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.score = None
        self.no_ip = None
        self.no_input = None
        self.seq = ''
    def to_bedstr(self):
        ''' convert to bedtool compatible string format'''
        bedstr = '{} {} {} {} {} {} {} {}\n'.format(self.chrom, self.start, self.end,0,'.',self.strand, '.', '.')
        return bedstr

def enrich_score(sites, input_sites, cluster):
    ''' calculate binomial enrichment score for Cluster object 
        add attribute  `score`, `no_ip`, `no_input` into object'''
    ip_in_cluster = 0
    input_in_cluster = 0
    # find how many read from IP/Input in cluster
    for s in sites:
        
        if s>= cluster.start and s<= cluster.end:
            ip_in_cluster += 1
    for i in input_sites:
        if i>= cluster.start and i<= cluster.end:
            input_in_cluster += 1
    
    # probability of success
    p = len(sites)/(len(sites)+len(input_sites))
    
    
    # binom
    total_cluster = ip_in_cluster + input_in_cluster
    
    bn = binom(total_cluster, p)
    prob = bn.cdf(ip_in_cluster)
    cluster.score = prob
    cluster.no_ip = ip_in_cluster
    cluster.no_input = input_in_cluster

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

    # convert to numpy array
    if combine:
        X = np.array(sites + input_sites).reshape(1, -1).T
    else:
        X = np.array(sites).reshape(1, -1).T

    # run DBSCAN
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(X)

    # convert to Cluster Object
    lbl=clustering.labels_[clustering.core_sample_indices_]
    pos=X.flatten()[clustering.core_sample_indices_]
    d=pd.DataFrame([lbl, pos]).T
    dmax = d.groupby(0).max().values.flatten().tolist()
    dmin = d.groupby(0).min().values.flatten().tolist()

    clusters = []
    for start, end in zip(dmin, dmax):
        # create cluster object
        c = Cluster(interval.chrom, start, end, interval.strand)
        
        # calculate enrichment score
        enrich_score(sites, input_sites, c)

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
            control_cluster.append(Cluster(interval.chrom, start, end, interval.strand))
            cindex += 1
            lindex += 1
    
    return control_cluster

