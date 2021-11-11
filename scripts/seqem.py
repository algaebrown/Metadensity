import pandas as pd
import math
from .sequence import *
from sklearn.linear_model import LinearRegression
import math
import numpy as np
from RBPamp.seed import PSAMSetBuilder
from collections import Counter
import numpy as np
from scipy.stats import entropy
import seaborn as sns

import matplotlib.pyplot as plt

################################# Z-score ##############################

### is moved to sequence

############################## Determine score cutoff #############################
def plot_fit(fit, z_sort, min_index=0, score = 'z-score', save = None):
    f= plt.figure()
    plt.plot(fit, label = 'fit')
    plt.plot(z_sort.values, label = 'data')
    plt.xlabel('Rank')
    plt.ylabel(score)
    plt.scatter(min_index, z_sort.values[min_index], label = 'cutoff')
    plt.legend()
    plt.title('Determine {} cutoff'.format(score))
    
    if save:
        f.savefig(save)
    else:
        plt.show()
    plt.close(f)

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
#################################### RBP AMP PSAM #####################################
def get_kmer_set(zscore, z_score_thres = 15):
    ''' from kmer zscore df, return kmer_set for PSAM_builder'''
    s = zscore.sort_values(ascending = False)
    s = s[s>z_score_thres]
    kmer_set = []
    for kmer in s.index:
        if 'N' not in kmer:
            kmer_set.append((kmer, s[kmer]))
    return kmer_set


def best_score_window(seq, psam, best_score):
    ''' best best window scoring for cluster'''
    windows = [seq[i:i+psam.n] for i in range(len(seq)) if i+ psam.n < len(seq)]
    score = [psam.score(w)/best_score for w in windows]
    
    # find best scoring window
    if len(score) == 0:
        # seq is shorter than psam
        print(seq, psam.consensus)
        best_window_score = 0
        best_window_start = 0
    else:
        best_window_score = max(score)
        best_window_start = score.index(best_window_score)
    
    return [best_window_score, best_window_start]
def crosslink_site_relative_to_window(csln, strand, start, end, best_window_start):
    ''' calculate crosslink position relative to best window start
    csln = list of crosslink in genome coords
    strand: +,-
    start = cluster start in genome coords
    end = cluster end in genome coords
    best_window_start: start in seq (flanking-cluster-flanking)
    seq_window: length of flanking for seq'''
    if type(csln) == float:
        print(csln)
    
    if strand == '-':
        cluster_start = end
        relative_csln = [cluster_start-pos for pos in csln]
    else:
        cluster_start = start
        relative_csln = [pos-cluster_start for pos in csln] # distance to 5'
    
    #relative to window start
    relative_to_window = [pos-best_window_start for pos in relative_csln]
    
    return relative_to_window
    
from multiprocessing import Pool
def calculate_cluster_score(df, psams, psam_thres = 0.8, n_pool = 8, timeout = 1000):
    ''' master function to find best window and crosslink position for each cluster'''
    use_seq = []
    for i, psam in enumerate(psams):
        psam= psam.shrink(thresh = psam_thres)
        use_seq.append(psam.consensus)
        best_score = psam.score(psam.consensus)
        
        
        # calculate score
        pool = Pool(int(n_pool))
        tasks = [
            (
            s,
            psam,
            best_score
                    
            )
             for s in df['seq']
             ] # pseudocount
        jobs = [pool.apply_async(best_score_window, task) for task in tasks]
        list_of_series = []
        for job in jobs:
            if job.get(timeout=timeout):
            
                # if return anything
                list_of_series.append(job.get(timeout=timeout))
        pool.close()
        
        
        s = pd.DataFrame(list_of_series, columns = ['psam_{}_score'.format(i), 'psam_{}_start'.format(i)], index = df.index)
        
        df = pd.concat([df, s], axis = 1) # append into the same df
        
        
        # get crosslink site
              
        
        df['csln_site_{}'.format(i)] = df.apply(lambda x: crosslink_site_relative_to_window(x['ip_site'], x['strand'], x['start'], x['end'], x['psam_{}_start'.format(i)]), axis=1)
        
    return df, use_seq
    


def csln_site_to_dist(site, pos_min = -15, pos_max = 15):
    # count # site in position
    c = Counter([s for s in site if s>=pos_min and s<= pos_max])

    # make to probability
    p = counter_to_prob(c)
    
    
    return p.iloc[0].values # np.array
    
def crosslink_distribution(all_csln_sites, normalize = True):
    ''' each cluster weighted equally '''
    prob_from_each_site = [csln_site_to_dist(site) for site in all_csln_sites]
    #print([p.shape for p in prob_from_each_site if len(p)!=30])
    mean = np.mean(np.stack(prob_from_each_site), axis = 0)
    
    
    return mean
        
    
    
def counter_to_prob(counter,pos_min = -15, pos_max = 15):
    ''' count the # crosslink from each cluster'''
    r = list(range(pos_min, pos_max+1))
    prob = pd.DataFrame(columns = r)
    
    prob = prob.append(counter, ignore_index = True)
    
    prob.fillna(0, inplace = True)
    
    prob = prob + 1
    prob = prob/np.sum(prob.values)
    return prob


def learn_crosslink_preference(df,
                                pos_min = -15, pos_max = 15):
    ''' learn crosslink preference for each motif from the highest pseudoscore clusters'''
    # find best exaplained PSAM
    psams_cols = [c for c in df.columns if 'psam' in c and 'score' in c]
    df['best_explained'] = df[psams_cols].idxmax(axis = 1)
    
    all_p = pd.DataFrame(columns = list(range(pos_min, pos_max+1)))
    sub_df = df.loc[df['INCLUDE']==True]
    for name, group in sub_df.groupby(by = 'best_explained'):
        group_number = name.split('_')[1]
        # count crosslink occurence from the top supprotee clusters
        
        
        
        #c = Counter(np.concatenate(top_cluster['csln_site_{}'.format(group_number)].values))
        #make to probability
        #prob = counter_to_prob(c)
        
        p = crosslink_distribution(group['csln_site_{}'.format(group_number)].values)
        all_p.loc[name] = p
        
    return all_p
        
        
    
    

def likelihood(csln, prob):
    
    csln_dist = pd.DataFrame(columns = prob.index)
    site_count = Counter(csln)
    
    log_L = 0
    for pos in site_count.keys():
        if pos in prob.index:
            log_L += np.log(prob[pos])*site_count[pos]
    
    
    return log_L/len(csln)
def best_likelihood(row, all_p=None):
    psam_key = row['best_explained']
    group_number = psam_key.split('_')[1]
    if psam_key in all_p.index:
        L = likelihood(row['csln_site_{}'.format(group_number)],all_p.loc[psam_key])
    else:
        # None of the INCLUDE is best represented by that psam, unlikely to be true
        L = -100
    
    return L
    
def whole_likelihood(df, all_p):
    ''' calculate the likelihood for the whole dataframe'''
    all_L = []
    
    
    df['log_L'] = df.apply(best_likelihood, all_p = all_p, axis = 1)
def update_cluster_boundry(row, use_seq = None, psam_window = 3):
    ''' update start with csln_site[], 
    update end with clsn_site{} + psam len (depends on psam score!)'''
    if row['INCLUDE'] == True:
        psam_group = int(row['best_explained'].split('_')[1])
        
        # update sequence, and extend 3nt
        start = int(row['psam_{}_start'.format(psam_group)])-psam_window # psam_start and add 3nt window
        if start < 0:
            start = 0
        end = int(row['psam_{}_start'.format(psam_group)]) + len(use_seq[psam_group]) + psam_window # psam_start + psam_len+ 3nt window
        if end > len(row['seq']):
            end = len(row['seq'])
        row['seq'] = row['seq'][start: end]

        # update cluster boundry
        if row['strand'] == '+':
            row['start'] = row['start'] + start
            row['end'] = row['start'] + end
        else:
            row['end'] = row['end'] - start
            row['start'] = row['start'] - end        
    return row
    
def split_cluster():
    pass
    # TODO
    # if cluster length > 30?
    
    # if included, extract GCAUG
    
    # split IP, input sites (+- 10 window)
    
    # other parts treat as seperate cluster (in total max is 3)
def update(df, use_seq, log_L_thres = -3, cluster_stamp = 0):
    ''' update INCLUDE and the cluster boundry'''
    # record when it is included
    df.loc[(df['log_L']>log_L_thres) & (df['INCLUDE'].isnull()), 'STAMP'] = cluster_stamp
    df.loc[df['log_L']>log_L_thres, 'INCLUDE'] = True
    
    # update cluster boundry: only the newly added ones
    df.loc[df['STAMP']==cluster_stamp].apply(update_cluster_boundry, use_seq = use_seq, axis = 1)
        
    # delete columns psam_{}_score; psam_{}_start; csln_site_{}
    to_drop = [c for c in df.columns if 'psam' in c] + [c for c in df.columns if 'csln' in c]
    df.drop(labels = to_drop, inplace = True, axis = 1)
    
    return df
def kmer():
    pass
    
    