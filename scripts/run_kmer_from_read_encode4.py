from metadensity.kmer_from_read import *
import pandas as pd
import pysam
from scipy.stats import pearsonr
import sys
import os
import numpy as np
import math
from sklearn.linear_model import LinearRegression


rbp_file = pd.read_pickle('/home/hsher/projects/ClipNet/ENCODE_stats/rbp_df.pickle')


def return_fobj(uid, rbp_file = rbp_file):
    '''return BedTool, Pysam objects'''
    raw_path = '/home/hsher/seqdata/eclip_raw/'
    
    row = rbp_file.loc[rbp_file['uid']==uid]
    if row.shape[0] == 0:
        print('No matching data')
    else:
        bam1 = row['bam_0'].values[0]
        bam2 = row['bam_1'].values[0]
        bam_in1 = row['bam_control_0'].values[0]
        bam_in2 = row['bam_control_1'].values[0]
    
    
    bam1_fobj = pysam.Samfile(bam1, 'rb')
    bam2_fobj = pysam.Samfile(bam2, 'rb')
    bam_input1_fobj = pysam.Samfile(bam_in1, 'rb')
    bam_input2_fobj = pysam.Samfile(bam_in2, 'rb')
    
    return bam1_fobj, bam2_fobj, bam_input1_fobj, bam_input2_fobj
def zscore_cutoff(zscore, diff_cutoff = 1, plot = False, score = 'z-score', save = None):
    ''' based on linear regression, figure out what z-score threshold is rapidly increasing'''
    z_sort = zscore.sort_values()
    
    # use first half
    length = math.ceil(len(z_sort)/2)
    X = np.reshape(np.arange(length),(-1, 1))
    y = z_sort.values[:length]

    reg = LinearRegression().fit(X, y)
    
    fit = reg.predict(np.array(np.reshape(np.arange(len(z_sort)),(-1, 1))))
    diff = z_sort.values - fit
    
    # positions where the real score is greater than background
    surpass_background = np.where(diff > diff_cutoff)[0]
    surpass_background = surpass_background[surpass_background > length] # must be greater than median
    
    if len(surpass_background) > 0:
        min_index = np.min(surpass_background)
        zcutoff = z_sort.values[min_index]
        if plot:
            plot_fit(fit, z_sort, min_index, score = score, save = save) 
    else:
        zcutoff = None
        # No one is significant
        plot_fit(fit, z_sort, save = save) 

    return zcutoff
def filtered_corr(z):
    ''' only significant z-score is considered'''
    z.dropna(inplace = True)
    thres1 = zscore_cutoff(z[0])
    thres2 = zscore_cutoff(z[1])
    
    fail1 = z[0][z[0]< thres1].index # kmers
    fail2 = z[1][z[1]< thres2].index # kmers
    
    masked_z = z.copy()
    masked_z.loc[fail1, 0] = 0
    masked_z.loc[fail2, 1] = 0
    
    return masked_z

if __name__=='__main__':
    uid = sys.argv[1]
    outdir = '/home/hsher/eclip_read_kmer'
    
    rep1, rep2, sminput1, sminput2 = return_fobj(uid)
    # calculate z-score
    z_score = main(rep1, rep2, sminput1, bam_input2=sminput2, chrom =  'chr1', start = 0, end = 18956422, n_iter = 100, n_sample = 1000, single_end = True)

    z_score.to_pickle(os.path.join(outdir, '{}.zscore'.format(uid)))
    print(z_score.idxmax(axis = 0))

    # filter and get pearson corr
    filtered_zscore = filtered_corr(z_score)
    filtered_zscore.dropna(inplace = True)
    filtered_zscore.replace(np.inf, filtered_zscore.max().max())
    r = pearsonr(filtered_zscore[0], filtered_zscore[1])[0]
    print('pearson corr: {:.2f}'.format(r))

    with open(os.path.join(outdir, 'pearson.csv'), 'a') as f:
        f.write('{},{}\n'.format(uid, r))