#'''
#Created on Jul 28, 2020
#@author: Hsuan-lin Her
#'''

import sys
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter
from .truncation import *
from .sequence import *
import matplotlib.pyplot as plt
from collections import Counter
from optparse import OptionParser
from multiprocessing import Pool
import pandas as pd
from pybedtools import BedTool
import pysam
import math
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
    def make_fixed_size(self, size = 30):
        ''' make cluster fixed sized '''
        if self.end-self.start > size:
            print('Warning: cluster size is {}'.format(self.end-self.start))
        
        center = math.ceil((self.start + self.end)/2)
        self.start = center - math.ceil(size/2)
        self.end = self.start + size

    def enrich_score(self, total_ip_in_region, total_input_in_region):
        ''' calculate enrichment score by binom'''
    
        
        # probability of success
        p = (total_ip_in_region+1)/(total_ip_in_region + total_input_in_region+2)
    
        
        # binom
        total_cluster = self.no_ip + self.no_input
    
        bn = binom(total_cluster, p)
        prob = bn.cdf(self.no_ip)
        self.score = prob
        
        # pseudocount
        
        bn = binom(total_cluster+2, p)
        prob = bn.cdf(self.no_ip+1)
        self.pseudoscore = prob
        
        print(total_ip_in_region, self.no_ip, total_input_in_region, self.no_input, prob)
    
    
    def to_bedstr(self):
        ''' convert to bedtool compatible string format'''
        bedstr = '{} {} {} {} {} {} {} {}\n'.format(self.chrom, self.start, self.end,0,'.',self.strand, '.', '.')
        return bedstr


def find_cluster(interval, bam_fileobj, inputbam_fileobj, combine = False, eps = 2, min_samples = 2, size = 30, read2 = True):
    ''' find read start cluster using DBSCAN
    interval: BedTool interval
    bam_fileobj: pysam IP bam
    inputbam_fileobj: pysam Input bam
    combine: if True, create cluster using both IP and Input reads; better for compare to Input.
    eps, min_samples: DBSCAN param
    '''
    # fetch read starts
    sites = read_start_sites(bam_fileobj, interval = interval, read2 = read2)
    input_sites = read_start_sites(inputbam_fileobj, interval = interval, read2 = read2)
    
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
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(X) #TODO epi and min_samples needs to be tuned

    # convert to Cluster Object
    
    # not all belong to a cluster
    lbl=clustering.labels_[clustering.core_sample_indices_] # labels
    pos=X.flatten()[clustering.core_sample_indices_] # positions
    identity = identity[clustering.core_sample_indices_] # IP or Input
    
    d=pd.DataFrame([lbl, pos, identity], index = ['label', 'position', 'identity']).T
    dmax = d.groupby('label')['position'].max().values.flatten().tolist()
    dmin = d.groupby('label')['position'].min().values.flatten().tolist()
    
    # TODO: break down gigantic cluster
    ip_members = [g[1].loc[g[1]['identity'] == 'IP', 'position'].tolist() for g in d.groupby('label')]
    input_members = [g[1].loc[g[1]['identity'] == 'Input', 'position'].tolist() for g in d.groupby('label')]

    clusters = []
    for start, end, ip_site, input_site in zip(dmin, dmax, ip_members, input_members):
        # create cluster object
        c = Cluster(interval.chrom, start, end, interval.strand, ip_site, input_site)
        c.make_fixed_size(size = size)
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



def get_cluster_seq(clusters):
    ''' find cluster sequence'''
    seqs = []
    for c in clusters:
        
        if type(c.start)!= int or type(c.end) != int:
            print(c.start, c.end)
        s = get_interval_seq(c.chrom, c.start, c.end, c.strand)
        c.seq = s
def region_cluster(interval, bam, inputbam, combine = True, eps = 2, min_samples = 2, size = 30, vis = False, read2 = True):
    
    ''' find cluster within bedtools interval '''
    
    # create file object
    bam_fileobj = pysam.Samfile(bam, 'rb')
    inputbam_fileobj = pysam.Samfile(inputbam, 'rb')
    
    # find clusters combining input and IP read start
    clusters, input_sites, sites = find_cluster(interval, bam_fileobj,inputbam_fileobj, combine = combine, eps = eps, min_samples = min_samples, size = size, read2 = read2)
    # write the sequence into cluster object
    get_cluster_seq(clusters)
    
    # generate size-matched, region-matched control cluster
    control_clusters = control_cluster(clusters, interval)
    get_cluster_seq(control_clusters)
    
    if vis:
        visualize_cluster(input_sites, sites, clusters, control_clusters)
    
    # convert result to list
    clus_list = [[c.chrom, c.start, c.end, c.strand, c.no_ip, c.no_input, c.pseudoscore, str(c.seq), c.ip_site, c.input_site]
                 for c in clusters]
    ctrl_list = [[c.chrom, c.start, c.end, c.strand, c.no_ip, c.no_input, c.score, str(c.seq)] 
                 for c in control_clusters]
    
    return clus_list, ctrl_list

def region_cluster_ri(interval, bam, inputbam, pseudocount = 0.1, sigma = 2, ri_height = 0.02, size = 30, read2 = True):
    ''' use RI to center the cluster instead of DBSCAN '''
    # create file object
    bam_fileobj = pysam.Samfile(bam, 'rb')
    inputbam_fileobj = pysam.Samfile(inputbam, 'rb')
    
    ip=np.array(truncation_relative_axis(bam_fileobj, interval = interval, read2 = read2))
    in_=np.array(truncation_relative_axis(inputbam_fileobj, interval = interval, read2 = read2))

    total_ip_in_region = ip.sum()
    total_input_in_region = in_.sum()
    
    
    ipdist = ip+pseudocount
    ipdist = gaussian_filter(ipdist/ipdist.sum(), sigma = sigma)
    
    indist = in_+pseudocount
    indist = gaussian_filter(indist/indist.sum(), sigma = sigma)
    trun = ipdist*np.log(ipdist/indist) ### relative entropy

    # get sequence
    seq = getRNAsequence(interval)
    # ip peaks
    peaks, _ = find_peaks(trun, height=ri_height, distance = size/2)
    
    def get_cluster_info(peak):
        ''' local function to get info for clusters '''
        clus_list = []
        for p in peaks:
            chrom = interval.chrom
            start = int(p-size/2)
            end = int(p+size/2)
            strand = interval.strand

            absl_start = interval.start + start if strand == '+' else interval.end-end
            absl_end = interval.start + end if strand == '+' else interval.end-start

            
            # how many reads
            no_ip = ip[start:end].sum()
            no_in = in_[start:end].sum()

            # probability of success
            p = (total_ip_in_region+1)/(total_ip_in_region + total_input_in_region+2)
        
            # binom
            total_cluster = no_ip + no_in
            bn = binom(total_cluster+2, p)
            pseudoscore = bn.cdf(no_ip+1)
            
            subseq = str(seq[start:end])
            
            clus_list.append([chrom, absl_start, absl_end, strand, no_ip, no_in, pseudoscore, subseq])
        return clus_list
    
    clus_list = get_cluster_info(peaks)
    
    peaks, _ = find_peaks(-trun, height=0, distance = size/2)
    ctrl_list = get_cluster_info(peaks)

    
    
    return clus_list, ctrl_list
    

def main(ip_bamfile, input_bamfile, enrich_bed = None, n_upper = None, n_pool = 8, timeout = 1000, eps = 2, min_samples = 2, combine = True, size = 30, tsv = None, use_ri =False, read2 = True):
    ''' DB scan for cluster
    ip_bamfile:
    input_bamfile:
    enrich_bed: bed file containing enriched region (produced by region_call.py)
    '''
    
    if enrich_bed:
        regions = BedTool(enrich_bed)
    else:
        # Evan's region caller
        df = pd.read_csv(tsv, sep = '\t')
        df = df.loc[(df['q_betabinom_1'] < 0.05) & (df['q_betabinom_2'] < 0.05)]
        regions = BedTool.from_dataframe(df[['chr', 'start', 'end', 'details', 'gc','strand']])
    # maximal regions
    if n_upper:
        regions = BedTool(regions[:n_upper]).saveas()
    
    print('Segmented Regions')

    # setup multiprocessing
    if not use_ri:
        pool = Pool(int(n_pool)) #interval, ip_bamfile, input_bamfile, ip_total, smi_total, pval_thres = 0.05
        tasks = [
                (
                interval,
                ip_bamfile,
                input_bamfile,
                combine,
                eps,
                min_samples,
                size, 
                read2
                
                )
                for interval in regions
                ] # pseudocount
        print('doing DB-SCAN for {} regions:'.format(len(tasks)))
        jobs = [pool.apply_async(region_cluster, task) for task in tasks]
    else:
        pseudocount = 0.1
        sigma = 2
        ri_height = 0.02
        size = 30 #TODO allow options for these
        pool = Pool(int(n_pool)) #interval, ip_bamfile, input_bamfile, ip_total, smi_total, pval_thres = 0.05
        tasks = [
                (
                interval,
                ip_bamfile,
                input_bamfile,
                pseudocount,
                sigma,
                ri_height,
                size,
                read2
                
                )
                for interval in regions
                ] # pseudocount
        print('doing DB-SCAN for {} regions:'.format(len(tasks)))
        jobs = [pool.apply_async(region_cluster_ri, task) for task in tasks]

    clusters = []
    control_clusters = []
    
    for job in jobs:
        if job.get(timeout=timeout):
            
            # if return anything
            clusters.extend(job.get(timeout=timeout)[0])
            control_clusters.extend(job.get(timeout=timeout)[1])
    
    # make dataframe
    print('Clusters found', len(clusters))
    print('Control clusters matched', len(control_clusters))
    
    try:
        clus_df = pd.DataFrame(clusters, 
                            columns = ['chrom', 'start', 'end', 'strand', 
                                        'no_ip', 'no_input', 'pseudoscore', 'seq', 'ip_site', 'input_site'])
    except:
        clus_df = pd.DataFrame(clusters, 
                            columns = ['chrom', 'start', 'end', 'strand', 
                                        'no_ip', 'no_input', 'pseudoscore', 'seq'])
    
    # control df
    ctrl_df = pd.DataFrame(control_clusters, 
                           columns = ['chrom', 'start', 'end', 'strand', 'no_ip', 'no_input', 'score', 'seq'])
    # should never overlap
    ctrl_bed = BedTool.from_dataframe(ctrl_df[['chrom', 'start','end', 'score', 'seq', 'strand', 'no_ip', 'no_input']])
    real_bed = BedTool.from_dataframe(clus_df[['chrom', 'start','end', 'pseudoscore', 'seq', 'strand']])

    valid_ctrl = ctrl_bed.intersect(real_bed, s = True, v = True).saveas()
    ctrl_df = valid_ctrl.to_dataframe(names = ['chrom', 'start','end', 'score', 'seq', 'strand', 'no_ip', 'no_input'])
    
    return clus_df, ctrl_df

