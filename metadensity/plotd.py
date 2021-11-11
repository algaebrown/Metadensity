''' Created on June, 2020
@author:Hsuan-lin Her
Scripts to visualize meta-density
'''
from collections import OrderedDict
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys
from .metadensity import *
from .pos_enrich import PosEnrichResult
import math
from . import settings
#################### generate axis ###################################
featnames = ['exon', 'intron']
ax_width_dict = settings.ax_width_dict
feat_len_dict = settings.feat_len

generic_rna = ['first_exon', 'exon','intron', 'last_exon']
protein_coding = ['five_prime_UTR', 'first_CDS', 'CDS', 'last_CDS', 'three_prime_UTR']
branchpoints = ['branchpoint', 'branchpoint_pred']
polyAs = ['polyAsignal', 'polyAsite']

def calcaulte_grid_width(features_to_show, ax_width_dict):
    ''' calculate the number of grid needed for gridspec '''
    width = {}
    
    exist_keys = list(ax_width_dict.keys())
    
    for feat in features_to_show:
        for key in exist_keys:
            if key in feat: # if string 'CDS' in 'first_CDS'
                valid_key = key
                break
        try:
            width[feat] = ax_width_dict[valid_key]
        except:
            width[feat] = ax_width_dict[feat]

    return width


def generate_axis(nrows = 2, ax_width_dict = ax_width_dict, color_bar = False, feat_len_dict = feat_len_dict, features_to_show = featnames, anno = None, height = None):
    ''' generate matplotlib ax for feature plots '''
    
    width_dict = calcaulte_grid_width(features_to_show, ax_width_dict)
    
    total_width = sum(list(width_dict.values()))*2
    if color_bar:
        total_width += 1
    
    if height is None:
        height = nrows*3
    
    
    
    # generate figure
    fig = plt.figure(figsize = (total_width, height))
    spec = gridspec.GridSpec(ncols=total_width, nrows=nrows, figure=fig)
    
    ax_dict = {}
    
    for r in range(nrows):
        current = 0
        for feat in features_to_show:
            for align in ['left', 'right']:
                # search for the key in ax_width_dict
                
                width = width_dict[feat]
                
                flen = feat_len_dict[feat]
                
                
                ## left y-axis
                if feat == features_to_show[0] and align == 'left':
                    ax_dict[feat, align, 'rep{}'.format(r+1)] = fig.add_subplot(spec[r,current: current+width])
                else:
                    ax_dict[feat, align, 'rep{}'.format(r+1)] = fig.add_subplot(spec[r,current: current+width], sharey = ax_dict[features_to_show[0], 'left', 'rep{}'.format(r+1)])
                    plt.setp(ax_dict[feat, align, 'rep{}'.format(r+1)].get_yticklabels(), visible=False)
                    
                if r == 0: # if first row, add title
                    if align == 'left':
                        lbl = "5\'"
                    else:
                        lbl = "3\'"

                    if feat == 'five_prime_UTR':
                        featstr = '5\' UTR'
                    elif feat == 'three_prime_UTR':
                        featstr = '3\' UTR'
                    else:
                        featstr=feat
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_title(f'{featstr}({lbl})')
                if r+1 < nrows: # if not last row, set xticklabels invisible
                    plt.setp(ax_dict[feat, align, 'rep{}'.format(r+1)].get_xticklabels(), visible=False)
                
                                  
                
                
                ## X-ticks with special feature
                
                if align == 'right':
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen+1, flen/5))
                    xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
                    
                    if feat == 'three_prime_UTR':
                        xticklabel[-1] = 'TTS'
                    if feat == 'intron':
                        xticklabel[-1] = '3\' SS'
                    if feat == 'last_CDS':
                        xticklabel[-1] = 'stop codon'
                    if 'branchpoint' in feat or 'polyA' in feat: # point features
                        ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen, flen/5))
                        xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
                        xticklabel[0] = feat


                    
                else:
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen, flen/5))
                    xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
                    if feat == 'five_prime_UTR':
                        xticklabel[0] = 'TSS'
                    if feat == 'intron':
                        xticklabel[0] = '5\' SS'
                    if feat == 'first_CDS':
                        xticklabel[0] = 'start codon'
                    if 'branchpoint' in feat or 'polyA' in feat: # point features
                        ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen+1, flen/5))
                        xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
                        xticklabel[-1] = feat
                        
                
                ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticklabels(xticklabel, rotation = 90)    
                current = current+width
                
    
    if color_bar:
        ax_dict['colorbar'] = fig.add_subplot(spec[:,-1:])
    
    return fig, ax_dict


    

################################## real plotting functions start here #############################################
def plot_rbp_map(metas, alpha = 0.6, ymax = 0.001, features_to_show = featnames, sort = False, rep_handle = 'combined'):
    ''' get a bunch of Metadensity or Metatruncation Object, plot their individual density in a heatmap'''
    fig, ax_dict = generate_axis(nrows = len(metas), color_bar = True, features_to_show = features_to_show)
    
    # set ylabel
    
    _ = [ax_dict[key].set_ylabel('transcripts') for key in ax_dict.keys() if 'five_utr' in key and 'left' in key]
    
   
    for m, rep in zip(metas, ['rep{}'.format(i) for i in range(1, len(metas)+1)]):
        i=0
        
        # autodetect object type
        if isinstance(m, Metatruncate):
            den_arr = m.density_array
        elif isinstance(m, Metadensity):
            den_arr = m.density_array
        else:
            print('Feeding wrong object {}. Only accept Metadensity or Metatruncate'.format(type(m)))
        
        # plotting features
        for feat in features_to_show:
            for align in ['left', 'right']:
                # extract density
                if rep_handle == 'combined':
                    # mean of all reps
                    density_concat = np.nanmean(np.stack([den_arr[feat,align, r] for r in m.eCLIP.rep_keys]), axis = 0)
                elif rep_handle == 'all':
                    # plot each rep individually
                    density_concat = np.concatenate([den_arr[feat,align, r] for r in m.eCLIP.rep_keys], axis = 0)
                else:
                    # specific rep
                    density_concat = den_arr[feat,align, rep_handle]
                
                # sort lexicographically for better visualization
                if sort:
                    density_concat = density_concat[np.lexsort(np.rot90(density_concat))]
                
                # find axis and plot
                ax = ax_dict[feat, align, rep]
                sns.heatmap(density_concat, ax = ax, cbar_ax = ax_dict['colorbar'], vmin = 0, vmax = ymax, cmap="YlGnBu")
                
                if align == 'left' and feat == features_to_show[0]:    
                    ax.set_ylabel(m.name)
                    
                i+= 1
    #plt.suptitle('RBP map: individual transcript')
    return fig



# show that std is large
def plot_mean_density(metas, ymax = 0.001, alpha = 0.3, plot_std = True, stat = 'mean', features_to_show = featnames, smooth = False, color_dict = None, mask = False):
    ''' get a bunch of eCLIPs, plot their mean density'''
    fig, ax_dict = generate_axis(nrows = 1,  features_to_show = features_to_show)

    
    #_ = [ax_dict[key].set_ylabel(stat + ' density') for key in ax_dict.keys() if features_to_show[0] in key and 'left' in key]
    _ = [ax_dict[key].set_ylim(ymax = ymax, ymin = 0) for key in ax_dict.keys() if features_to_show[0] in key]
    
    for m in metas:
        i=0
        if isinstance(m, Metatruncate):
            den_arr = m.density_array
        elif isinstance(m, Metadensity):
            den_arr = m.density_array
        elif isinstance(m, PosEnrichResult):
            if mask:
                den_arr = m.masked
            else:
                den_arr = m.stat
        else:
            print('Feeding wrong object {}. Only accept Metadensity or Metatruncate'.format(type(m)))

        for feat in features_to_show:
            for align in ['left', 'right']:

                ####### make replicate  together and mean #########
                density_concat = np.concatenate([den_arr[feat,align, r] for r in m.eCLIP.rep_keys], axis = 0)
                flen = feat_len_dict[feat]
                if align == 'left':
                    density_concat = density_concat[:, :flen]
                else:
                    density_concat = density_concat[:, -flen:]

                ####### calculating mean and median #########
                if stat == 'mean':
                    md = np.nanmean(density_concat, axis = 0)
                if stat == 'median':
                    md = np.nanmedian(density_concat, axis = 0)
                
                std = np.nanstd(density_concat, axis = 0)
                n = density_concat.shape[0]
                sem = std/np.sqrt(n)

                ####### get the right axis #########
                ax = ax_dict[feat, align, 'rep1']
                
                

                if smooth:
                    if color_dict:
                        ax.plot(gaussian_smooth(md), label = m.name, color = color_dict[m.name])
                    else:
                        ax.plot(gaussian_smooth(md), label = m.name)
                else:
                    if color_dict:
                        ax.plot(md, label = m.name, color = color_dict[m.name])
                    else:
                        ax.plot(md, label = m.name)
                    
                    if plot_std:
                        if color_dict:
                            ax.fill_between(np.arange(len(md)), md-sem, md+sem, label = m.name, alpha = alpha, color = color_dict[m.name])
                        else:
                            ax.fill_between(np.arange(len(md)), md-sem, md+sem, label = m.name, alpha = alpha)
                
                if align == 'left' and feat == features_to_show[0]:
                    if m.background_method == 'relative information':
                        ylbl = '{} relative information'.format(stat)
                    elif m.background_method == None:
                        ylbl = '{} IP'.format(stat)
                    else:
                        ylbl = '{} subtracted '.format(stat)
                    
                    if m.normalize:
                        ylbl += ' normalized density'
                    else:
                        ylbl += ' raw'
                    ax.set_ylabel(ylbl)
                
                
                    
                i+= 1
    ncol = math.ceil(len(metas)/10)
    plt.legend(bbox_to_anchor = (1.5, 1), ncol = ncol)
    
    
    #plt.suptitle('Pooled Density from Many Transcript')
    return fig

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def auto_rbp_color(rbp_list):
    cmap = plt.cm.get_cmap('rainbow', len(rbp_list)+1)
    return dict(zip(rbp_list, [cmap(i) for i in range(len(rbp_list))]))

def plot_enrich(metas, ymax = 0.5, alpha = 0.3, features_to_show = featnames, smooth = False, color_dict = None, mask = False, scatter = False, ymin =0):
    ''' get a bunch of eCLIPs, plot their mean density'''
    fig, ax_dict = generate_axis(nrows = 1,  features_to_show = features_to_show)

    #_ = [ax_dict[key].set_ylabel(stat + ' density') for key in ax_dict.keys() if features_to_show[0] in key and 'left' in key]
    _ = [ax_dict[key].set_ylim(ymax = ymax, ymin = ymin) for key in ax_dict.keys() if features_to_show[0] in key]
    
    if not color_dict:
        rbp_list = [m.name for m in metas]
        # automatically generate colors
        color_dict = auto_rbp_color(rbp_list)

    for m in metas:
        if isinstance(m, PosEnrichResult):
            if mask:
                den_arr = m.masked
            else:
                den_arr = m.test_stat
        else:
            print('Feeding wrong object {}. Only accept PosEnrich'.format(type(m)))

        for feat in features_to_show:
            for align in ['left', 'right']:
                ####### get the right axis #########
                ax = ax_dict[feat, align, 'rep1']

                for rep in m.rep_keys:
                    md = den_arr[feat, align, rep]
                    if smooth:
                        if scatter:
                            ax.scatter(list(range(len(md))), md, label = m.name, color = color_dict[m.name], alpha = alpha, marker = '.')
                    
                        ax.plot(gaussian_smooth(md), label = m.name, color = color_dict[m.name])
                        
                    else:
                        
                        ax.plot(md, label = m.name, color = color_dict[m.name])
                        
                
                if align == 'left' and feat == features_to_show[0]:
                    
                    ylbl = '{} statistic'.format(m.test)
                    
                    ax.set_ylabel(ylbl)
    # remove redundant legends
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ncol = math.ceil(len(metas)/10)
    plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor = (1.5, 1), ncol = ncol)
    #plt.suptitle('Positional Enrichment Result')
    return fig
