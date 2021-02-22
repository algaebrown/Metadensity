import numpy as np
from scipy.stats import ks_2samp
from .truncation import  truncation_relative_axis, read_start_sites
from .metadensity import gaussian_smooth
import deepdish as dd

def read_icshape(fname):
    data = {}
    with open(fname) as f:
        for line in f:
            
            values = line[:-1].split('\t') # remove \n
            gene_id = values[0]
            gene_len = values[1]
            coverage = values[2]
            reactivity = [float(v) if v!= 'NULL' else np.nan for v in values[3:]]
            data[gene_id] = reactivity
    return data

def window_around(data, start, end):
    ''' return values from start:end, handles problems like start < 0 or end > data langth '''
    values = [np.nan] * (end-start)
    if start < 0:
        # need to pad zero in front
        value_starts_at = (-start)
        start = 0
    else:
        value_starts_at = 0
        
    if end > len(data):
        value_ends_at = len(data) - end
        end = len(data)
    else:
        value_ends_at = len(values)
    
    d = data[start:end]
    values[value_starts_at:value_ends_at] = d
    
    return values
    

def window_around_eclip_sites(interval, data, bam, window = 10, single_end = False):
    ''' for a bedtool interval, return windows of shape values'''
    eclip = truncation_relative_axis(bam, interval = interval, single_end = single_end)[:-1]
    shape = np.array(data[interval.attrs['ID']])
    
    all_windows = []
    for i in range(len(eclip)):
        if eclip[i]>0: # when there is truncation
            
                
            window_val = window_around(shape, i-window, i+window)
            if np.sum(np.isnan(window_val)) < len(window_val):
                for n in range(eclip[i]): 
                    all_windows.append(window_val)
    if len(all_windows) > 0:
        return np.stack(all_windows)
    else:
        return None
            
def many_intervals(bam, data, filtered, num = 50, window = 10, single_end = False):
    '''for many genes that have shape data, return all the windows'''
    all_w = []
    for interval in filtered[:num]:
        w = window_around_eclip_sites(interval, data, bam, window = window, single_end = single_end)
        if w is not None:
            all_w.append(w)
    return np.concatenate(all_w, axis = 0)

def windows_for_all(bam1, bam2, bam_in1, data, filtered, bam_in2 = None, window = 10,  num = 50, single_end = False):
    w_ip1 = many_intervals(bam1, data, filtered, window = window, num = num, single_end = single_end)
    w_ip2 = many_intervals(bam2, data, filtered, window = window, num = num, single_end = single_end)
    w_in =  many_intervals(bam_in1, data, filtered, window = window, num = num, single_end = single_end)
    
    if bam_in2:
        w_in2 = many_intervals(bam_in2, data, filtered, window = window, num = num, single_end = single_end)
        return w_ip1, w_ip2, w_in, w_in2
    else:
        return w_ip1, w_ip2, w_in

def ks_all_pos(w, w_in, p_thres = 0.01):
    ''' run KS test for all position '''
    ps = []
    kss = []
    for pos in range(w.shape[1]):
        
        # if ip is larger than in
        ks_large, pval_large = ks_2samp(w[:, pos], w_in[:, pos], alternative = 'less')
        
        # if ip is smaller than in
        ks_small, pval_small = ks_2samp(w[:, pos], w_in[:, pos], alternative = 'greater')
        
        if pval_large < pval_small:
            ps.append(pval_large)
            kss.append(ks_large)
            
        else:
            ps.append(pval_small)
            kss.append(-ks_small)
            
    #masked = np.ma.masked_where(np.array(ps)>p_thres, np.array(kss)) masked array cannot be saved
    return ps, kss
