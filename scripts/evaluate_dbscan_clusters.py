import sys
import os
from metadensity.sequence import *
from metadensity.truncation import *
from pybedtools import BedTool
from dataloader import *
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp, entropy
def df_to_bed(df):
    ''' preprocess new peak'''
    df.iloc[:, :8].drop_duplicates(inplace = True) # the "list" is not usable
    df['fold'] = (df['no_ip']+1)/(df['no_input']+1)
    
        
    try:
        # for peaks
        bedcol= ['chrom', 'start', 'end', 'no_ip', 'no_input', 'strand', 'pseudoscore', 'fold']
        return BedTool.from_dataframe(df[bedcol])
    except:
        # for control regions
        bedcol =  ['chrom', 'start', 'end', 'no_ip', 'no_input', 'strand', 'score', 'fold']
        return BedTool.from_dataframe(df[bedcol])

def purity_sensitivity(new_peak, old_peak):
    ''' given dbscan output, bedtool overlap with IDR peaks '''
    purity = len(new_peak.intersect(old_peak, s = True, u = True))/len(new_peak)
    sensitivity = len(old_peak.intersect(new_peak, s = True, u = True))/len(old_peak)

    return purity, sensitivity
    
# fetch crosslinking distribution for all kmers
def get_kmer_pos(motif, start, end, strand, seq):
    ''' given seq, chrom. strat end strand, return the starting point of a given motif'''
    try:
        start_from_five_prime = seq.index(motif)
    except:
        return None, None
        
    if strand == '-':
        
        end_motif = end-start_from_five_prime
        start_motif = end_motif - len(motif)
        
    else:
        start_motif = start + start_from_five_prime
        end_motif = start_motif + len(motif)
    return start_motif, end_motif

def crosslink_pos(df, motif, bam_ip, bam_in, window = 30, single_end = False):
    ''' return truncation count in IP and IN library for every site in the df '''
    ins = []
    ips = []
    for index, row in df.iterrows():
        motif_start, motif_end = get_kmer_pos(motif, row['start'], row['end'], row['strand'], row['seq'])
        if motif_start is not None: # if motif exist in seq
           
            trun_ip = truncation_relative_axis(bam_ip, chrom=row['chrom'], start=motif_start-window, end=motif_end+window, strand=row['strand'], single_end = single_end)
            trun_in = truncation_relative_axis(bam_in, chrom=row['chrom'], start=motif_start-window, end=motif_end+window, strand=row['strand'], single_end = single_end)
        
            ins.append(trun_in)
            ips.append(trun_ip)
    return ins, ips

def crosslink_entropy(ins, ips, stat = 'median'):
    ''' calculate crosskinking relative entropy '''
    
    ip_csln_den = (np.stack(ips)+1)/(np.sum(np.stack(ips)+1, axis = 1))[:, None]
    in_csln_den = (np.stack(ins)+1)/(np.sum(np.stack(ins)+1, axis = 1))[:, None]

    # mean is driven by extreme values
    if stat == 'mean':
        ent = entropy(ip_csln_den.mean(axis = 0), in_csln_den.mean(axis = 0)) 
    # median is also strongly affected by lots of crosslink in IN
    elif stat == 'median':
        med_ip = np.median(ip_csln_den, axis = 0)
        med_in = np.median(in_csln_den, axis = 0)
        ent = entropy(med_ip/med_ip.sum(), med_in/med_in.sum())
    return ent

def motif_ks(df, motif):
    ''' calculate motif pseudoscore KS p-value '''
    df['motif'] = df['seq'].str.contains(motif)

    contain_motif = df.loc[df['motif']]
    no_motif = df.loc[df['motif'] == False]
    return ks_2samp(contain_motif['pseudoscore'], no_motif['pseudoscore'], alternative = 'less') # test for enriched
    
def KS_pseudoscore_all_kmers(df, z, score = 'AUC'):
    ''' tryout high zscore kmers'''
    
    
    scores_ = pd.Series(index = z.index)
    
    for motif in z.index:
        ks,pval = motif_ks(df, motif)
            
        scores_[motif] = pval
    return scores_


def get_kmer_stat(df, df_control, bam_ip, bam_in, outfile, k = 6, pval_thres = 0.001, window = 30, n_sample = 500, single_end = False, pseudoscore_thres = 0.7):
    ''' the main function '''
    ######### k-mer zscore ###########

    mean, std = simulate_kmer_background(df_control['seq'].tolist(), k = k, n_sample = df.loc[df['pseudoscore']>pseudoscore_thres,'seq'].shape[0])
    z = kmer_zscore(df.loc[df['pseudoscore']>pseudoscore_thres,'seq'], mean, std, k = k) # need to filter for high scoring clusters


    ######### KS test for pseudoscore ######
    print('calculating pseudocode enrichment for top kmers')
    ks_pvals = KS_pseudoscore_all_kmers(df, z.sort_values()[-top_kmer:], score = 'KS')


    ######### to file ######
    print('combineind data')
    kmer_stats = pd.concat([z, ks_pvals], axis = 1)
    kmer_stats.columns = ['zscore', 'pseudoscore_ks_pval']

    print('top pseudoscore kmers: ')
    print(kmer_stats.sort_values(by = 'pseudoscore_ks_pval').iloc[:5])

    ######### get crosslinking bias, relative entropy ######
    # find significant k-mers
    signi = kmer_stats.loc[kmer_stats['pseudoscore_ks_pval']<pval_thres]
    if signi.shape[0] > 10:
        signi = signi.sort_values(by = 'pseudoscore_ks_pval').iloc[:10]
    for kmer in signi.index:
        with_kmer = df.loc[df['seq'].str.contains(kmer)]
        if with_kmer.shape[0]>n_sample:
            with_kmer = with_kmer.sample(n_sample)
        ins, ips = crosslink_pos(with_kmer, kmer, bam_ip, bam_in, window = window, single_end = single_end) ### TODO, maybe consider bam2 as well
        relative_entropy = crosslink_entropy(ins, ips)
        kmer_stats.loc[kmer, 'csln_entropy'] = relative_entropy
    
    if 'csln_entropy' in kmer_stats.columns:
        # some may not have any significant kmers
        print('top direct kmers: ')
        print(kmer_stats.sort_values(by = 'csln_entropy', ascending = False).iloc[:5])

    kmer_stats.to_csv(outfile)

if __name__ == '__main__':
    print(sys.argv[1])
    df = pd.read_pickle(sys.argv[1])
    df_control = pd.read_pickle(sys.argv[2])
    outdir =sys.argv[3]
    print(df.shape, df_control.shape)

    prefix = os.path.basename(sys.argv[1]).split('.')[0]
    rep = os.path.basename(sys.argv[1]).split('.')[1]
    uid = prefix.split('_')[0] # 204

    ######## load bams ############
    try:
        # ENCODE 3 data
        bam1, bam2, bamin = return_fobj3(uid)
        single_end = False
        encode4 = False
    except:
        # ENCODE 4 data
        bam1, bam2, bamin, bamin2 = return_fobj4(uid)
        single_end = True
        encode4 = True
        if rep == 'rep2':
            bamin=bamin2

    ####### Calculate k-mer #######
    for k in [3,4,5,6,7,8]:
        top_kmer = 100
        outfile = os.path.join(outdir, prefix+'.{}.{}mer.csv'.format(rep, k))

        if rep=='rep1':
            get_kmer_stat(df, df_control, bam1, bamin, outfile, k = k, single_end = single_end)
        else:
            get_kmer_stat(df, df_control, bam2, bamin, outfile, k = k, single_end = single_end)


        ######## Compare the top clusters against the worst clusters ##########
        outfile = os.path.join(outdir, prefix+'.{}.{}mer.topvsworse.csv'.format(rep, k))
        if rep=='rep1':
            get_kmer_stat(df.loc[df['pseudoscore']>0.6], df.loc[df['pseudoscore']<0.2], bam1, bamin, outfile, k = k, single_end = single_end)
        else:
            get_kmer_stat(df.loc[df['pseudoscore']>0.6], df.loc[df['pseudoscore']<0.2], bam2, bamin, outfile, k = k, single_end = single_end)

    ######## Compare against CLIPper outputs ##########
    if encode4:
        idr = BedTool(encode4_data.loc[encode4_data['uid']==uid, 'idr'].iloc[0])
        if rep == 'rep1':
            indv = BedTool(encode4_data.loc[encode4_data['uid']==uid, 'bed_0'].iloc[0])
        else:
            indv = BedTool(encode4_data.loc[encode4_data['uid']==uid, 'bed_1'].iloc[0])
    else:
        idr = BedTool(encode_data.loc[encode_data['uid']==uid, 'idr'].iloc[0])
        if rep == 'rep1':
            indv = BedTool(encode_data.loc[encode_data['uid']==uid, 'bed_0'].iloc[0])
        else:
            indv = BedTool(encode_data.loc[encode_data['uid']==uid, 'bed_1'].iloc[0])
    # real clusters
    df['pseudoscore_bin'] = pd.cut(df['pseudoscore'], bins = [0,0.2,0.4, 0.6, 0.8, 1])
    stats = []
    for name, group in df.groupby(by = 'pseudoscore_bin'):
        size = group.shape[0]
        pur, sen = purity_sensitivity(df_to_bed(group), idr)
        pur_indv, sen_indv = purity_sensitivity(df_to_bed(group), indv)
        stats.append([name, pur, sen, pur_indv, sen_indv, size])
    
    # control clusters
    pur, sen = purity_sensitivity(df_to_bed(df_control), idr)
    pur_indv, sen_indv = purity_sensitivity(df_to_bed(df_control), indv)
    size = df_control.shape[0]
    stats.append(['ctrl', pur, sen, pur_indv, sen_indv, size])
    
    stat_df = pd.DataFrame(stats, columns = ['pseudoscore_group','purity_idr', 'sen_idr','purity_indv', 'sen_indv', 'n_clus'])
    outfile = os.path.join(outdir, prefix+'{}.stat.csv'.format(rep))
    stat_df.to_csv(outfile)

    print(stat_df)

    

