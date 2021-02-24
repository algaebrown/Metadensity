from .metadensity import *
from .plotd import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

def highly_exp_biogps(cell_line = 'HEPG2', transcript_type = 'protein_coding', sample_no = 200):
    ''' pre-construct metagene for highly expressed transcripts '''
    biogps = pd.read_csv('/home/hsher/projects/Metadensity/notebooks/gene_attribute_matrix.txt', sep = '\t')
    biogps = biogps.iloc[2:]
    biogps.iloc[:, 2:] = biogps.iloc[:, 2:].astype(float, copy = False)
    
    highexp_genes = biogps.loc[biogps[cell_line]>0, '#'].tolist()
    tids = [t.attrs['transcript_id'] for t in transcript if t.attrs['transcript_type'] == transcript_type and t.attrs['gene_name'] in highexp_genes] # hungtingtin
    cds_metagenes = Build_many_metagene(tids, sample_no = sample_no)

    return cds_metagenes

def construct_distribution(e, metagenes, use_truncate = True):
    ''' Given eCLIP object, construct null and IP density '''

    ############# NULL ####################
    if use_truncate:
        m_null=Metatruncate(e, e.name+'_null', metagenes = metagenes, background_method = 'get null', normalize = True, smooth = False)
        m_null.get_density_array(use_truncation = True)
        print('Done constructing an empirical null distribution from Input')

        ############# IP ######################
        m_ip=Metatruncate(e, e.name+'_real', metagenes = metagenes, background_method = None, normalize = True, smooth = False)
        m_ip.get_density_array(use_truncation = True)
        print('Done constructing IP distribution')
    else:
        m_null=Metadensity(e, e.name+'_null', metagenes = metagenes, background_method = 'get null', normalize = True, smooth = False)
        m_null.get_density_array()
        print('Done constructing an empirical null distribution from Input')

        ############# IP ######################
        m_ip=Metadensity(e, e.name+'_real', metagenes = metagenes, background_method = None, normalize = True, smooth = False)
        m_ip.get_density_array()
        print('Done constructing IP distribution')

    return m_null, m_ip
############################# wilcox test method ##########################################
from scipy.stats import wilcoxon
def auc_cum_wilcox(real, null, n_largest_to_remove = 5):
    ''' return AUC for a feature throughout all positions'''
    
    
    all_diff = []
    all_p = []
    all_ks = []
    for pos in range(real.shape[1]):
        
        # sort and minus
        pos_real = np.sort(real[:, pos])
        pos_null = np.sort(null[:, pos])
        
        #diff = pos_real - pos_null
        #all_diff.append(np.mean(diff[:-n_largest_to_remove])) # mean difference, remove largest few values prevent outlier

        # TODO mask KS insignificant
        if np.sum(pos_real[:-n_largest_to_remove] - pos_null[:-n_largest_to_remove]) == 0:
            all_p.append(1)
            all_ks.append(0)
        else:
            ks, pval = wilcoxon(pos_real[:-n_largest_to_remove], pos_null[:-n_largest_to_remove], alternative = 'greater')
            all_p.append(pval)
            all_ks.append(ks)
    return np.array(all_ks), np.array(all_p)
def Wilcox_enrich(meta_null, meta_ip, pval_thres = 0.05, n_largest_to_remove = 5):
    ''' main function to calculate enrichment based on AUC? '''
    result = {}
    pvals = {}
    for key in meta_null.truncate_array.keys():
        real = meta_ip.truncate_array[key]
        real = np.nan_to_num(real, 0)
        null = meta_null.truncate_array[key]
        null = np.nan_to_num(null, 0)

        auc_per_pos, all_p = auc_cum_wilcox(real, null, n_largest_to_remove = n_largest_to_remove)

        # TODO multiple hypothesis testing correction
        #result[key] = np.ma.masked_where(all_p < 0.05, auc_per_pos)
        result[key] = auc_per_pos
        pvals[key] = all_p
    return result, pvals

############################# KS-test method ##########################################
def auc_cum_dist(real, null, n_largest_to_remove = 5):
    ''' return AUC for a feature throughout all positions'''
    
    
    all_diff = []
    all_p = []
    all_ks = []
    for pos in range(real.shape[1]):
        
        # sort and minus
        pos_real = np.sort(real[:, pos])
        pos_null = np.sort(null[:, pos])
        
        #diff = pos_real - pos_null
        #all_diff.append(np.mean(diff[:-n_largest_to_remove])) # mean difference, remove largest few values prevent outlier

        # TODO mask KS insignificant
        ks, pval = ks_2samp(pos_real[:-n_largest_to_remove], pos_null[:-n_largest_to_remove], alternative = 'lesser')
        all_p.append(pval)
        all_ks.append(ks)
    return np.array(all_ks), np.array(all_p)

def AUC_enrich(meta_null, meta_ip, pval_thres = 0.05, n_largest_to_remove = 5):
    ''' main function to calculate enrichment based on AUC? '''
    result = {}
    pvals = {}
    for key in meta_null.truncate_array.keys():
        real = meta_ip.truncate_array[key]
        real = np.nan_to_num(real, 0)
        null = meta_null.truncate_array[key]
        null = np.nan_to_num(null, 0)

        auc_per_pos, all_p = auc_cum_dist(real, null, n_largest_to_remove = n_largest_to_remove)

        # TODO multiple hypothesis testing correction
        #result[key] = np.ma.masked_where(all_p < 0.05, auc_per_pos)
        result[key] = auc_per_pos
        pvals[key] = all_p
    return result, pvals

########################## percentile reject #####################################

def null_thres(null, perc = 99):
    ''' return threshold of each position that need to be rejected'''
     
    # calculate 99th percentile threshold
    thres = np.percentile(null, perc, axis = 0)

    # thres cannot equal 0, interpolate using the closest non-zero value
    thres[thres==0]=np.nan
    ts = pd.Series(thres)
    ts=ts.interpolate(method='nearest').ffill().bfill()

    return ts.values

def mask_by_null(real, thres):
    ''' mask the density if not rejected'''
    
    # get stacked null distribution
    n_sample = real.shape[0]
    stacked_null= np.stack([list(thres)]*n_sample)
    
    # masking
    masked_real = np.ma.masked_where(real < stacked_null, real)
    np.ma.set_fill_value(masked_real, 0)
    return masked_real

def enrich_by_thres(meta_null, meta_ip):
    ''' main function for percentile threshold enrichment '''
    result = {}
    for key in meta_null.truncate_array.keys():
        real = meta_ip.truncate_array[key]
        real = np.nan_to_num(real, 0)
        null = meta_null.truncate_array[key]
        null = np.nan_to_num(null, 0)

        # get thres
        thres = null_thres(null)

        # get masked array
        masked_real = mask_by_null(real, thres)
        result[key] = masked_real
        
    return result

