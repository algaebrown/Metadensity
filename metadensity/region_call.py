'''
Created on Jul 28, 2020
@author: Hsuan-lin Her

'''


from optparse import OptionParser
import pandas as pd
from pybedtools import BedTool
from scipy.stats import poisson
from multiprocessing import Pool
import pysam
import numpy as np
def poission_pval(smi_total, smi_region, ip_total, ip_region):
    rate = smi_region / smi_total
    x = ip_total
    mu = rate*x
    
    
    # get pvalue
    return 1-poisson.cdf(ip_region, mu = mu)

def nreads(bam, interval):
    ''' pysam fetch but strand specific, return number reads mapped to region'''
    reads = bam.fetch(interval.chrom, interval.start, interval.stop)
    if interval.strand == '-':
        nreads = len([read for read in reads if read.is_reverse])
    else:
        nreads = len([read for read in reads if not read.is_reverse])
    return nreads

def region_based_calling(interval, ip_bamfile, input_bamfile, ip_total, smi_total, pval_thres = 0.05):
    ''' determine whether feature interval i is enriched in IP '''
    
    ip_bam = pysam.AlignmentFile(ip_bamfile)
    input_bam = pysam.AlignmentFile(input_bamfile)
    
    
    ip_region = nreads(ip_bam, interval)
    smi_region = nreads(input_bam, interval)

    # pseudocount: cause a lot of noise
    #ip_region += 1
    #smi_region += 1
    #ip_total += 1
    #smi_total += 1
    
    if ip_region > 0 and smi_region > 0:
        # if more than 1 reads in IP, check reads in SMInput
        
        # if more than 1 reads in SMInput
        pval = poission_pval(smi_total, smi_region, ip_total, ip_region)
        fold_enrich = np.log2((ip_region/ip_total)/(smi_region/smi_total))
        #print(ip_region, smi_region, ip_total, smi_total, pval, fold_enrich)

        if pval < pval_thres:
            # if unadjusted pvalue is significant, return
            try:
                bed = [interval.chrom, interval.start,interval.stop,interval.attrs['ID'],".",interval.strand, pval, fold_enrich, ip_region, smi_region]
            except:
                bed = [interval.chrom, interval.start,interval.stop,interval[3],".",interval.strand, pval, fold_enrich, ip_region, smi_region]
            return bed

def main(ip_bamfile, input_bamfile, features, n_upper = None, n_pool = 8, timeout = 1000, pval = 10**(-3), fold_change = 3, save_pickle = False):
    ''' run_poission given eCLIP object'''
    
    # fetch total reads
    ip_bam = pysam.AlignmentFile(ip_bamfile)
    ip_total = ip_bam.mapped
    smi_bam = pysam.AlignmentFile(input_bamfile)
    smi_total = smi_bam.mapped
    print('Total reads in IP: {}; Input {} '.format(ip_total, smi_total))

      
    feat = BedTool(features)
    # maximal regions
    if n_upper:
        feat = BedTool(feat[:n_upper]).saveas()
    else:
        n_upper = len(feat)
    print('done extract features')

    # setup multiprocessing
    pool = Pool(int(n_pool)) #interval, ip_bamfile, input_bamfile, ip_total, smi_total, pval_thres = 0.05
    tasks = [
            (
            interval,
            ip_bamfile,
            input_bamfile,
            ip_total,
            smi_total
            
            )
             for interval in feat
             ] # pseudocount
    print('calling for {} regions:'.format(len(tasks)))

    unadjusted = []
    jobs = [pool.apply_async(region_based_calling, task) for task in tasks]
    for job in jobs:
        if job.get(timeout=timeout):
            
            # if return anything
            unadjusted.append(job.get(timeout=timeout))
    
    # make dataframe
    print('no feature unadj p-val < 0.05', len(unadjusted))
    df = pd.DataFrame(unadjusted, columns = ['chrom', 'start', 'stop', 'ID', '.', 'strand', 'pval', 'log2_fold', 'ip_region', 'smi_region'])
    
    # bonferrni correction

    # find second smallest pval to replace 0, avoid log10(0)
    second_smallest_pval = df['pval'].drop_duplicates().nsmallest(2).iloc[-1]
    df.loc[df['pval'] == 0, 'pval'] = second_smallest_pval
    print('replace pval = 0 with 2nd smallest {}'.format(second_smallest_pval))
    df['log10_padj'] = np.log10(df['pval']* n_upper)
    
    # filter
    filtered = df.loc[(df['log10_padj']< np.log10(pval)) & (df['log2_fold'] > np.log2(fold_change))]
    print('after filter:', filtered.shape)

    # BedTool header 
    bedcol= ['chrom', 'start', 'stop', 'ID', '.', 'strand', 'log10_padj', 'log2_fold']
    
    if save_pickle:
        return BedTool.from_dataframe(filtered[bedcol]), df
    else:
        return BedTool.from_dataframe(filtered[bedcol])

def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS CLIPPER FOR REGIONAL ENRICHMENT
        python region_call.py --ip <bamfile> --input <bamfile> --feature <gff> --out <fname>
        """
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012, Hsuan-lin Her 2020.
                         CLIP peakfinder by regional enrichment."""
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--ip", "-i", dest="IPbam", help="IP bam file", type="string", metavar="FILE.bam")

    parser.add_option("--input", "-s", dest="SMInputbam", help="size match input bam file")
    
    parser.add_option("--feature", "-f", dest="feature", help="gff3 file containing regions")
    parser.add_option("--out", "-o", dest="outfile", help="output file path")
    parser.add_option("--maxfeat", "-n", dest="nfeat", help="maximal number of feature", type = "int")
    parser.add_option("--thread", "-t", dest="thread", help="multiprocessing", default = 8, type = "int")
    parser.add_option("--pval", "-p", dest="pval", help="p-value threshold", default = 0.001, type = "float")
    parser.add_option("--fold", "-r", dest="fold", help="fold change threshold", default = 8, type = "float")
    parser.add_option("--pickle",  dest="pickle", help="save original dataframe", action = "store_true")
    
    return parser

if __name__ == "__main__":
    parser = option_parser()
    (options, args) = parser.parse_args()
    # call main
    if options.pickle:
        enriched_bed, df = main(options.IPbam, options.SMInputbam, options.feature, n_upper = options.nfeat, 
        n_pool = options.thread, timeout = 1000, pval = options.pval, fold_change = options.fold, save_pickle = options.pickle)
        df.to_pickle(options.outfile.replace('bed', 'pickle'))
    else:
        enriched_bed = main(options.IPbam, options.SMInputbam, options.feature, n_upper = options.nfeat, 
        n_pool = options.thread, timeout = 1000, pval = options.pval, fold_change = options.fold)
    # save file
    enriched_bed.saveas(options.outfile)

