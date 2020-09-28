'''
Created on Jul 28, 2020
@author: Hsuan-lin Her
'''

import sys

from truncation import *
from sequence import *
import matplotlib.pyplot as plt
from collections import Counter
from optparse import OptionParser
from multiprocessing import Pool
import pandas as pd
from pybedtools import BedTool
import pysam

def get_cluster_seq(clusters, window = 10):
    ''' find cluster sequence with window'''
    seqs = []
    for c in clusters:
        
        if type(c.start)!= int or type(c.end) != int:
            print(c.start, c.end)
        s = get_interval_seq(c.chrom, c.start-window,c.end+window, c.strand)
        c.seq = s
def region_cluster(interval, bam, inputbam, vis = False, combine = True, window = 10):
    ''' find cluster within bedtools interval '''
    
    # create file object
    bam_fileobj = pysam.Samfile(bam, 'rb')
    inputbam_fileobj = pysam.Samfile(inputbam, 'rb')
    
    # find clusters combining input and IP read start
    clusters, input_sites, sites = find_cluster(interval, bam_fileobj,inputbam_fileobj, combine = combine)
    # write the sequence into cluster object
    get_cluster_seq(clusters, window = window)
    
    # generate size-matched, region-matched control cluster
    control_clusters = control_cluster(clusters, interval)
    get_cluster_seq(control_clusters, window = window)
    
    if vis:
        visualize_cluster(input_sites, sites, clusters, control_clusters)
    
    # convert result to list
    clus_list = [[c.chrom, c.start, c.end, c.strand, c.no_ip, c.no_input, c.pseudoscore, str(c.seq), c.ip_site, c.input_site]
                 for c in clusters]
    ctrl_list = [[c.chrom, c.start, c.end, c.strand, c.no_ip, c.no_input, c.score, str(c.seq)] 
                 for c in control_clusters]
    
    return clus_list, ctrl_list

def main(ip_bamfile, input_bamfile, enrich_bed, n_upper = None, n_pool = 8, timeout = 1000):
    ''' DB scan for cluster
    ip_bamfile:
    input_bamfile:
    enrich_bed: bed file containing enriched region (produced by region_call.py)
    '''
    
         
    regions = BedTool(enrich_bed)
    # maximal regions
    if n_upper:
        regions = BedTool(regions[:n_upper]).saveas()
    
    print('Segmented Regions')

    # setup multiprocessing
    pool = Pool(int(n_pool)) #interval, ip_bamfile, input_bamfile, ip_total, smi_total, pval_thres = 0.05
    tasks = [
            (
            interval,
            ip_bamfile,
            input_bamfile
            
            
            )
             for interval in regions
             ] # pseudocount
    print('doing DB-SCAN for {} regions:'.format(len(tasks)))

    clusters = []
    control_clusters = []
    jobs = [pool.apply_async(region_cluster, task) for task in tasks]
    for job in jobs:
        if job.get(timeout=timeout):
            
            # if return anything
            clusters.extend(job.get(timeout=timeout)[0])
            control_clusters.extend(job.get(timeout=timeout)[1])
    
    # make dataframe
    print('Clusters found', len(clusters))
    print('Control clusters matched', len(control_clusters))
    
    clus_df = pd.DataFrame(clusters, 
                           columns = ['chrom', 'start', 'end', 'strand', 
                                      'no_ip', 'no_input', 'pseudoscore', 'seq', 'ip_site', 'input_site'])
    ctrl_df = pd.DataFrame(control_clusters, 
                           columns = ['chrom', 'start', 'end', 'strand', 'no_ip', 'no_input', 'score', 'seq'])
    
    return clus_df, ctrl_df

def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS CLIPPER FOR REGIONAL SHRINKING
        python shrink_region.py --ip <bamfile> --input <bamfile> --bed <bed> --out <fname>
        
        To shrink enriched regions into binding sites that we hope to contain motifs
Returns clusters (putative binding sites) and control clusters (size-matched, region-matched).
By comparing cluster to control clusters, we can calculated R values.
Clustering algorithm is DBSCAN
        """
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012, Hsuan-lin Her 2020.
                         CLIP peakfinder by shrinking enriched region."""
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--ip", "-i", dest="IPbam", help="IP bam file", type="string", metavar="FILE.bam")

    parser.add_option("--input", "-s", dest="SMInputbam", help="size match input bam file")
    
    parser.add_option("--bed", "-b", dest="enrichBed", help="bedfile output by region_call.py")
    parser.add_option("--out", "-o", dest="outfile", help="output file prefix: /home/hsher/rbfox2")
    parser.add_option("--maxfeat", "-n", dest="nfeat", help="maximal number of feature", type = "int")
    parser.add_option("--thread", "-t", dest="thread", help="multiprocessing", default = 8, type = "int")
    
    
    return parser

if __name__ == "__main__":
    parser = option_parser()
    (options, args) = parser.parse_args()
    # call main
    
    
    clus_df, ctrl_df = main(options.IPbam, options.SMInputbam, options.enrichBed, n_upper = options.nfeat, 
    n_pool = options.thread, timeout = 1000)
    
    # save file
    clus_df.to_pickle(options.outfile + '.pickle')
    ctrl_df.to_pickle(options.outfile + '.ctrl.pickle')
