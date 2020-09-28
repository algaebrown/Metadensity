''' Created on June, 2020
@author:Hsuan-lin Her
Scripts to visualize meta-density
'''

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/hsher/projects/Metadensity')
from metadensity.metadensity import *

#################### generate axis ###################################

ax_width_dict = {'utr':1,
                'exon':2,
                'intron':3,
                'CDS':2}
feat_len_dict = {'five_utr':100, 
                 'three_utr':150, 
                 'intron': 1500, 
                 'exon':150, 
                 'CDS':150}
def calcaulte_grid_width(features_to_show, ax_width_dict):
    ''' calculate the number of grid needed for gridspec '''
    n_exon = len([f for f in features_to_show if 'exon' in f])
    n_cds = len([f for f in features_to_show if 'CDS' in f])
    n_utr = len([f for f in features_to_show if 'utr' in f])
    n_intron = len([f for f in features_to_show if 'intron' in f])

    width = (ax_width_dict['utr']*n_utr + ax_width_dict['exon']*n_exon + ax_width_dict['intron']*n_intron + ax_width_dict['CDS']*n_cds)*2
    return width


def generate_axis(nrows = 2, ax_width_dict = ax_width_dict, color_bar = False, feat_len_dict = feat_len_dict, features_to_show = featnames):
    ''' generate matplotlib ax for feature plots '''
    
    total_width = calcaulte_grid_width(features_to_show, ax_width_dict)
    if color_bar:
        total_width += 1
    height = nrows*3
    
    
    
    # generate figure
    fig = plt.figure(figsize = (total_width, height))
    spec = gridspec.GridSpec(ncols=total_width, nrows=nrows, figure=fig)
    
    ax_dict = {}
    
    for r in range(nrows):
        current = 0
        for feat in features_to_show:
            for align in ['left', 'right']:
                width = ax_width_dict[feat.split('_')[-1]]
                flen = feat_len_dict[feat.replace('first_','').replace('last_','')]
                
                
                ## left y-axis
                if feat == features_to_show[0] and align == 'left':
                    ax_dict[feat, align, 'rep{}'.format(r+1)] = fig.add_subplot(spec[r,current: current+width])
                else:
                    ax_dict[feat, align, 'rep{}'.format(r+1)] = fig.add_subplot(spec[r,current: current+width], sharey = ax_dict[features_to_show[0], 'left', 'rep{}'.format(r+1)])
                    plt.setp(ax_dict[feat, align, 'rep{}'.format(r+1)].get_yticklabels(), visible=False)
                    
                if r == 0: # if first row, add title
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_title(feat + '-'+ align)
                if r+1 < nrows: # if not last row, set xticklabels invisible
                    plt.setp(ax_dict[feat, align, 'rep{}'.format(r+1)].get_xticklabels(), visible=False)
                
                                  
                
                
                ## X-ticks with special feature
                
                if align == 'right':
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen+1, flen/5))
                    xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
                    
                    if feat == 'three_utr':
                        xticklabel[-1] = 'TTS'
                    if feat == 'intron':
                        xticklabel[-1] = '3\' SS'
                    if feat == 'last_CDS':
                        xticklabel[-1] = 'stop codon'
                    
                else:
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen, flen/5))
                    xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
                    if feat == 'five_utr':
                        xticklabel[0] = 'TSS'
                    if feat == 'intron':
                        xticklabel[0] = '5\' SS'
                    if feat == 'first_CDS':
                        xticklabel[0] = 'start codon'
                
                ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticklabels(xticklabel, rotation = 90)    
                current = current+width
                
    
    if color_bar:
        ax_dict['colorbar'] = fig.add_subplot(spec[:,-1:])
    return fig, ax_dict

####################### smoothing #######################################
def gaussian_smooth(mean, sigma = 5):
    ''' use gaussian kernel to smooth spiky data, sigma = stf '''
    y_vals = mean
    x_vals = np.arange(len(mean))
    
    smoothed_vals = np.zeros(y_vals.shape)
    for x_position in x_vals:
        kernel = np.exp(-(x_vals - x_position) ** 2 / (2 * sigma ** 2))
        kernel = kernel / sum(kernel)
        smoothed_vals[x_position] = sum(y_vals * kernel)
    return smoothed_vals
    


def plot_rbp_map(metas, quantile = False,  scaled = False, truncation = False, alpha = 0.6, ymax = 0.003, features_to_show = featnames):
    ''' get a bunch of eCLIPs, plot their individual density in a heatmap'''
    fig, ax_dict = generate_axis(nrows = len(metas), color_bar = True, features_to_show = features_to_show)
    
    # set ylabel
    
    _ = [ax_dict[key].set_ylabel('transcripts') for key in ax_dict.keys() if 'five_utr' in key and 'left' in key]
    
   
    for m, rep in zip(metas, ['rep{}'.format(i) for i in range(1, len(metas)+1)]):
        i=0
        if quantile:
            den_arr = m.qdensity_array
        if truncation:
            den_arr = m.truncate_array
        else:
            den_arr = m.density_array
        if scaled:
            den_arr = m.scaled_density_array
        for feat in features_to_show:
            for align in ['left', 'right']:
                # extract density
                density_concat = np.concatenate([den_arr[feat,align, r] for r in m.eCLIP.rep_keys], axis = 0)
                # sort lexicographically for better visualization
                #density_concat = density_concat[np.lexsort(np.rot90(density_concat))]
                
                ax = ax_dict[feat, align, rep]
                
                
                sns.heatmap(density_concat, ax = ax, cbar_ax = ax_dict['colorbar'], vmin = 0, vmax = ymax, cmap="YlGnBu")
                
                if align == 'left' and feat == 'five_utr':
                    
                    ax.set_ylabel(m.name)
                    
                i+= 1
    
    
    if quantile:
        t = 'Quantile Metadensity'
    elif scaled:
        t = 'Scaled Metadensity'
    elif truncation:
        t = 'truncation density'
    else:
        t = 'Normalized Metadensity'
    
    plt.suptitle(t)



def plot_density(eCLIP, 
                 five_utr_len=100, three_utr_len=150, intron_len = 1500, exon_len = 150, 
                 ymax = 0.03, logy = False, example = 'positive', quantile = False,scaled = False,
                alpha = 0.5):
    ''' plot eCLIP.density_array or eCLIP.qdensity_array individual level'''
    
    fig, ax_dict = generate_axis(nrows = 2)
    
    # density
    if quantile:
        eCLIP_den = eCLIP.qdensity_array
        title = 'Quantile Density'
    if scaled:
        eCLIP_den = eCLIP.scaled_density_array
        title = 'Scaled Density'
    else:
        eCLIP_den = eCLIP.density_array
        title = 'Normalized Density'
    
    
    # set xlim for exon, intron, utr's maximal length
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = three_utr_len) for key in ax_dict.keys() if 'three_utr' in key]
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = five_utr_len) for key in ax_dict.keys() if 'five_utr' in key]
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = exon_len) for key in ax_dict.keys() if 'exon' in key]
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = intron_len) for key in ax_dict.keys() if 'intron' in key]
    
    
    # set ylim for peak height
    _ = [ax_dict[key].set_ylim(ymin = 0, ymax = ymax) for key in ax_dict.keys() if 'five_utr' in key]
    _ = [ax_dict[key].set_ylabel(title) for key in ax_dict.keys() if 'five_utr' in key and 'left' in key]
    if logy:
        _ = [ax_dict[key].set_yscale('log') for key in ax_dict.keys() if 'five_utr' in key]
    
    
    for key in ax_dict.keys():
        eCLIP_key = tuple([example] + list(key))
        ax_dict[key].plot(eCLIP_den[eCLIP_key].T, alpha = alpha, color = 'indianred')
        
        
    plt.suptitle('{} {} : {} examples'.format(title, eCLIP.name,example))

def plot_prob(eCLIP, eCLIP_prob, example):
    '''plot prob distribution'''
    fig, ax_dict = generate_axis(nrows = 2, color_bar = True)
    
    # cm
    cm = sns.cubehelix_palette(50, hue=0.05, rot=0, light=0.9, dark=0, as_cmap = True)
    vmax = np.log(1)
    vmin = np.log(0.001)
    
    
    
    for keys in eCLIP_prob.keys():
        if example in keys:
            ax = ax_dict[keys[1:]]
            prob = eCLIP_prob[keys]
            if 'three_utr' in keys and 'rep2' in keys:
                sns.heatmap(np.log(prob[1:, :]), cmap=cm, ax = ax, cbar_ax = ax_dict['colorbar'], vmin = vmin, vmax = vmax) 
            else:
                sns.heatmap(np.log(prob[1:, :]), cmap=cm, ax = ax, cbar = False, vmin = vmin, vmax = vmax) 
                bins = prob.shape[0]
    
    # set label\n",
    ax_dict['colorbar'].set_ylabel('log probability')
    _ = [ax_dict[keys].set_ylabel('density bins') for keys in ax_dict.keys() if 'five_utr' in keys and 'left' in keys]
   
    _ = [ax_dict[keys].invert_yaxis() for keys in ax_dict.keys() if 'five_utr' in keys and 'left' in keys]      
    _ = [ax_dict[keys].set_yticklabels(np.arange(1,bins)) for keys in ax_dict.keys() if 'five_utr' in keys]
        
        
        
    
    
    plt.suptitle('{}: {} example'.format(eCLIP.name, example))

# show that std is large
def plot_mean_density(metas, ymax = 0.003, quantile = False, truncation = False, scaled = False, alpha = 0.3, plot_std = True, stat = 'mean', features_to_show = featnames, smooth = False):
    ''' get a bunch of eCLIPs, plot their mean density'''
    fig, ax_dict = generate_axis(nrows = 1,  features_to_show = features_to_show)
    
    
    
    # set ylabel
    
    _ = [ax_dict[key].set_ylabel(stat + ' density') for key in ax_dict.keys() if 'five_utr' in key and 'left' in key]
    _ = [ax_dict[key].set_ylim(ymax = ymax, ymin = 0) for key in ax_dict.keys() if 'five_utr' in key]
    
    for m in metas:
        i=0
        if quantile:
            den_arr = m.qdensity_array
        elif truncation:
            den_arr = m.truncate_array
        else:
            den_arr = m.density_array
        if scaled:
            den_arr = m.scaled_density_array
        for feat in features_to_show:
            for align in ['left', 'right']:
                density_concat = np.concatenate([den_arr[feat,align, r] for r in m.eCLIP.rep_keys], axis = 0)
                if stat == 'mean':
                    md = np.nanmean(density_concat, axis = 0)
                if stat == 'median':
                    md = np.nanmedian(density_concat, axis = 0)
                
                std = np.nanstd(density_concat, axis = 0)
                n = density_concat.shape[0]
                sem = std/np.sqrt(n)
                
                ax = ax_dict[feat, align, 'rep1']
                
                

                if smooth:
                    ax.plot(gaussian_smooth(md), label = m.name)
                else:
                    ax.plot(md, label = m.name)
                    if plot_std:
                        ax.fill_between(np.arange(len(md)), md-sem, md+sem, label = m.name, alpha = alpha)
                
                
                    
                i+= 1
    plt.legend(bbox_to_anchor = (1.5, 0.5))
    
    if quantile:
        t = 'Quantile Metadensity'
    elif truncation:
        t = 'Truncation Density'
    elif scaled:
        t = 'Scaled Metadensity'
    else:
        t = 'Normalized Metadensity'
    title = t
    plt.suptitle(title)
    

def plot_one_density(eCLIP, transcript_index,
                 five_utr_len=100, three_utr_len=150, intron_len = 1500, exon_len = 150, 
                 ymax = 0.03, logy = False, example = 'positive', quantile = False,
                alpha = 0.5):
    ''' plot eCLIP.density_array or eCLIP.qdensity_array for only 1 transcript'''
    
    fig, ax_dict = generate_axis(nrows = 2)
    
    # set xlim for exon, intron, utr's maximal length
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = three_utr_len) for key in ax_dict.keys() if 'three_utr' in key]
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = five_utr_len) for key in ax_dict.keys() if 'five_utr' in key]
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = exon_len) for key in ax_dict.keys() if 'exon' in key]
    _ = [ax_dict[key].set_xlim(xmin = 0, xmax = intron_len) for key in ax_dict.keys() if 'intron' in key]
    
    
    # set ylim for peak height
    _ = [ax_dict[key].set_ylim(ymin = 0, ymax = ymax) for key in ax_dict.keys() if 'five_utr' in key]
    _ = [ax_dict[key].set_ylabel('Normalized density') for key in ax_dict.keys() if 'five_utr' in key and 'left' in key]
    if logy:
        _ = [ax_dict[key].set_yscale('log') for key in ax_dict.keys() if 'five_utr' in key]
    
    # density
    if quantile:
        eCLIP_den = eCLIP.qdensity_array
    else:
        eCLIP_den = eCLIP.density_array
    for key in ax_dict.keys():
        eCLIP_key = tuple([example] + list(key))
        ax_dict[key].plot(eCLIP_den[eCLIP_key].T[transcript_index,:], alpha = alpha, color = 'indianred')
        
        
    plt.suptitle('Normalized density {} : {} examples'.format(eCLIP.name,example))

################################# Plots for consistency between probability, mean or median #############################
## consistency by scatter plot
colors = dict(zip(featnames, ['maroon','royalblue', 'navy', 'seagreen', 'lightskyblue', 'tomato', 'darkviolet', 'indigo', 'mediumpurple']))

def prob_consistency(eclip_prob, name, quantile = False):
    ''' plot consistency for eclip_prob '''
    
    if quantile:
        ncol = 1
        examples = ['positive']
    else:
        ncol = 2
        examples = ['positive','negative']
        
    f,ax = plt.subplots(1,ncol)
    
    if quantile:
        ax = [ax]
    
    i = 0
    for e in examples:
        for f in ['five_utr', 'exon', 'intron', 'three_utr']:
            for a in ['left', 'right']:
                ax[i].scatter(eclip_prob[e, f, a, 'rep1'].ravel(), eclip_prob[e, f, a, 'rep2'].ravel(), color = colors[f], label = f, alpha = 0.6, marker = '+')
        
        ax[i].set_xlabel('rep1 prob')
        ax[i].set_ylabel('rep2 prob')
        ax[i].set_title(e)
        i+=1
    plt.suptitle('probability distribution consistency: {}'.format(name))
    plt.tight_layout()
    plt.legend(loc = 'right', bbox_to_anchor = [1.5, 0.5])
    #plt.yscale('log')
    #plt.xscale('log')

def mean_density(eCLIP, use_quantile = False, use_scaled = False):
    ''' get mean and median, stderr density '''
    mean = {}
    median = {}
    stderr = {}
    
    if use_quantile:
        den_array = eCLIP.qdensity_array
    elif use_scaled:
        den_array = eCLIP.scaled_density_array
    else:
        den_array = eCLIP.density_array
    
    for key in den_array.keys():
        median[key] = np.nanmedian(den_array[key], axis = 0)
        mean[key] = np.nanmean(den_array[key], axis = 0)
        n = den_array[key].shape[0]
        stderr[key] = np.nanstd(den_array[key], axis = 0)/np.sqrt(n)
    return mean, median, stderr

def mean_med_consistency(meta, use_quantile = False, use_scaled = False, ymax = 0.002):
    ''' plot the consistency of mean and median for density array'''
    
    nrow = 1   
    
    f, ax = plt.subplots(nrow,2, figsize = (6,nrow*3) , sharex = True, sharey = True)
    ax = ax.flatten()
    mean, median, stderr = mean_density(meta, use_quantile = use_quantile, use_scaled = use_scaled)
    i = 0
    
    for f in featnames:
        for a in ['left', 'right']:
                
            ax[i].scatter(mean[f, a, 'rep1'], mean[f, a, 'm.eCLIP.rep_keys()'], label = f, color = colors[f], marker = '+', alpha = 0.6)
            ax[i].set_title('mean')
            ax[i].set_xlim((0,ymax))
            ax[i].set_ylim((0,ymax))
            ax[i].set_xlabel('rep1')
            ax[i].set_ylabel('rep2')
                
                
            ax[i+1].scatter(median[f, a, 'rep1'], mean[f, a, 'rep2'], label = f ,color = colors[f], marker = '+', alpha = 0.6)
            ax[i+1].set_title('median')
            ax[i+1].set_xlim((0,ymax))
            ax[i+1].set_ylim((0,ymax))
            ax[i+1].set_xlabel('rep1')
        
    plt.legend(bbox_to_anchor = (2, 1))
    plt.suptitle('Reproducibility {}'.format(meta.name))
    
    
    
########################################### entropy based #######################################################
def plot_entropy(eCLIPs, eCLIP_probs):
    ''' get a bunch of eCLIPs, plot their mean density'''
    fig, ax_dict = generate_axis(nrows = 1)
    # set ylabel
    _ = [ax_dict[key].set_ylabel('relative entropy') for key in ax_dict.keys() if 'five_utr' in key and 'left' in key]
    
    
    for eCLIP, prob in zip(eCLIPs, eCLIP_probs):
        
        for feat in featnames:
            for align in ['left', 'right']:
                pos = np.mean(np.array([prob['positive', feat,align, r] for r in ['rep1', 'rep2']]), axis = 0)
                neg = np.mean(np.array([prob['negative', feat,align, r] for r in ['rep1', 'rep2']]), axis = 0)
                entro = density_array_entropy(pos, neg)
                ax_dict[feat, align, 'rep1'].plot(entro, label = eCLIP.name, alpha = 0.5)
                
                
    plt.legend()
    plt.suptitle('relative entropy: positive v.s. negative')
    
                
    
    