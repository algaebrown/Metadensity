''' Created on June, 2020
@author:Hsuan-lin Her
Scripts to visualize meta-density
'''
from collections import OrderedDict
from scipy.ndimage import gaussian_filter
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys
from metadensity.metadensity import *
from metadensity.pos_enrich import PosEnrichResult
from metadensity.config import settings
import math
# from . import settings
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
                xticks, xticklabel = get_xticklabels(feat, align, n_tick = 5)
                ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(xticks)
                ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticklabels(xticklabel, rotation = 90)   


                current = current+width
                
    
    if color_bar:
        ax_dict['colorbar'] = fig.add_subplot(spec[:,-1:])
    
    return fig, ax_dict

def get_xticklabels(feat, align, n_tick = 5):
    '''create xticklabels for each kind of features'''
    flen = feat_len_dict[feat]

    if align == 'right':
        xticks = np.arange(0, flen+1, flen/n_tick)
        xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
        
        if feat == 'three_prime_UTR':
            xticklabel[-1] = 'TTS'
        if 'intron' in feat.split('_')[-1]:
            xticklabel[-1] = '3\' SS'
        if feat == 'last_CDS':
            xticklabel[-1] = 'stop codon'
        if 'branchpoint' in feat or 'polyA' in feat: # point features
            xticks = np.arange(0, flen, flen/5)
            xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
            xticklabel[0] = feat
    else:
        xticks = np.arange(0, flen, flen/5)
        xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
        if feat == 'five_prime_UTR':
            xticklabel[0] = 'TSS'
        if 'intron' in feat.split('_')[-1]:
            xticklabel[0] = '5\' SS'
        if feat == 'first_CDS':
            xticklabel[0] = 'start codon'
        if 'branchpoint' in feat or 'polyA' in feat: # point features
            xticks = np.arange(0, flen+1, flen/5)
            xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
            xticklabel[-1] = feat
    return xticks, xticklabel


def beautify(fig, offset = 10, trim = True, left = True):
    ''' wrapper around sns.despine to make ugly python figures better-looking'''
    for i,ax in enumerate(fig.get_axes()):
        if i>0:
            sns.despine(offset=offset, trim=trim, left = left, ax = ax)
        else:
            sns.despine(offset=offset, trim=trim, ax = ax)
        plt.setp(ax.get_xticklabels(), rotation=90)
    return fig

def rep_merge_options(den_arr, feat, align, metaden_object, rep_handle):
    '''
    return the appropriate array depending how user wants to merge the replicate
    rep_handle: 'mean' to plot the mean of 2 reps. 'concat' to show 2 reps individually. specify rep keys like 'rep1', 'rep2' to show only that rep.
    '''
    if rep_handle == 'mean':
        # mean of all reps
        density_concat = np.nanmean(np.stack([den_arr[feat,align, r] for r in metaden_object.eCLIP.rep_keys]), axis = 0)
    elif rep_handle == 'concat':
        # plot each rep individually
        density_concat = np.concatenate([den_arr[feat,align, r] for r in metaden_object.eCLIP.rep_keys], axis = 0)
    else:
        # specific rep
        density_concat = den_arr[feat,align, rep_handle]
    return density_concat

################################## real plotting functions start here #############################################

def plot_rbp_map(metas, alpha = 0.6, ymax = 0.001, features_to_show = featnames, sort = False, rep_handle = 'mean', cmap = 'Greys'):
    ''' get a bunch of Metadensity or Metatruncation Object, plot their individual density in a heatmap
    metas: list of Metadensity or Metatruncate object
    alpha: transparency in plt.plot()
    ymax: the max value in plt.set_ylim()
    features_to_show: list of genomic/transcriptomic features to show. options include all feature names. You can also use pre-set combinations such as `generic_rna`, `protein_coding` etc. Use metadensity.density_array to see what is available
    sort: whether to sort RNAs (rows in RBPmap) "lexicographically". setting true will make RNA with similar binding pattern cluster together on a map
    rep_handle: 'mean' to plot the mean of 2 reps. 'concat' to show 2 reps individually. specify rep keys like 'rep1', 'rep2' to show only that rep.
    cmap: color map to use 
    '''
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
                density_concat = rep_merge_options(den_arr, feat, align, m, rep_handle)
                
                # sort lexicographically for better visualization
                if sort:
                    density_concat = density_concat[np.lexsort(np.rot90(density_concat))]
                
                # find axis and plot
                ax = ax_dict[feat, align, rep]
                sns.heatmap(density_concat, ax = ax, cbar_ax = ax_dict['colorbar'], vmin = 0, vmax = ymax, cmap=cmap)
                
                if align == 'left' and feat == features_to_show[0]:    
                    ax.set_ylabel(m.name)
                
                ## X-ticks with special featurel sns.heatmap overrids xticks
                xticks, xticklabel = get_xticklabels(feat, align, n_tick = 5)
                
                ax_dict[feat, align, rep].set_xticks(xticks)
                ax_dict[feat, align, rep].set_xticklabels(xticklabel, rotation = 90)
                    
                i+= 1
    #plt.suptitle('RBP map: individual transcript')
    return fig



# show that std is large
def plot_mean_density(metas, ymax = 0.001, alpha = 0.3, plot_std = True, stat = 'mean', 
features_to_show = featnames, smooth = False, color_dict = None, mask = False, ylabel = None,
sigma = 5, rep_handle = 'concat'):
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
                density_concat = rep_merge_options(den_arr, feat, align, m, rep_handle)
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
                        ax.plot(gaussian_filter(md, sigma = sigma), label = m.name, color = color_dict[m.name])
                    else:
                        ax.plot(gaussian_filter(md, sigma = sigma), label = m.name)
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
                    if ylabel:
                        ax.set_label(ylbl)
                    else:
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

def plot_enrich(metas, ymax = 0.5, alpha = 0.3, features_to_show = featnames, smooth = False, color_dict = None, mask = False, scatter = False, ymin =0, sigma = 5):
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
                    
                        ax.plot(gaussian_filter(md, sigma = sigma), label = m.name, color = color_dict[m.name])
                        
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

# show that std is large
from scipy.ndimage import gaussian_filter, gaussian_filter1d
def make_concat_xtick(all_meta, features_to_show = featnames):
    cur_pos = 0
    xticks=[]
    xticklabels=[]
    
    for feat in features_to_show:
        for align in ['left', 'right']:
            flen = all_meta[0].density_array[feat, align, 'rep1'].shape[1]
        

                
            if align == 'right':
                xtick = np.arange(0, flen+1, flen/5)
                xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
                
                

                if feat == 'three_prime_UTR':
                    xticklabel[-1] = 'TTS'
                if 'intron' in feat:
                    xticklabel[-1] = '3\' SS'
                if feat == 'last_CDS':
                    xticklabel[-1] = 'stop codon'
                if 'branchpoint' in feat or 'polyA' in feat: # point features
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen, flen/5))
                    xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
                    xticklabel[0] = feat



            else:
                xtick=np.arange(0, flen, flen/5)
                xticklabel = ['{:.0f}'.format(x) for x in np.arange(0,flen, flen/5)]
                
                if feat == 'five_prime_UTR':
                    xticklabel[0] = 'TSS'
                if 'intron' in feat:
                    xticklabel[0] = '5\' SS'
                if feat == 'first_CDS':
                    xticklabel[0] = 'start codon'
                if 'branchpoint' in feat or 'polyA' in feat: # point features
                    ax_dict[feat, align, 'rep{}'.format(r+1)].set_xticks(np.arange(0, flen+1, flen/5))
                    xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen/5)]
                    xticklabel[-1] = feat


            xticks.append(xtick+cur_pos)
            xticklabels += xticklabel
            
            cur_pos += flen
                
    
    
    return np.concatenate(xticks), xticklabels


def plot_mean_density_concat(metas, ymax = 0.001, alpha = 0.3, plot_std = True, stat = 'mean', 
features_to_show = featnames, smooth = False, color_dict = None, mask = False, ylabel = None,
                            sigma = 5, figsize = (7,3)):
    ''' get a bunch of eCLIPs, plot their mean density'''
    f,ax = plt.subplots(figsize = figsize)
    ax.set_ylim(0, ymax)
    
    
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
        
        # get values
        density_concat = np.concatenate([m.concat_density_array(features = features_to_show, rep = r) for r in m.eCLIP.rep_keys], axis = 0)
        
        if smooth:
            density_concat = np.apply_along_axis(arr=density_concat, axis = 1, func1d=gaussian_filter1d, **{'sigma': sigma})
        
        if stat == 'mean':
            md = np.nanmean(density_concat,axis = 0)
        if stat == 'median':
            md = np.nanmedian(density_concat,axis = 0)
        
        
        std = np.nanstd(density_concat, axis = 0)
        n = density_concat.shape[0]
        sem = std/np.sqrt(n)
        
        if color_dict:
            ax.plot(md, color = color_dict[m.name],label = m.name)
        else:
            ax.plot(md, label = m.name)
    
        if plot_std:
            if color_dict:
                ax.fill_between(np.arange(len(md)), md-sem, md+sem, label = m.name, alpha = alpha, color = color_dict[m.name])
            else:
                ax.fill_between(np.arange(len(md)), md-sem, md+sem, label = m.name, alpha = alpha)
        
        
        
        
        
        if ylabel:
            ax.set_label(ylbl)
        else:
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

                
    # xticks
    xticks, xticklabels = make_concat_xtick(metas, features_to_show = features_to_show)
    print()
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation = 90)
                
    
    ax.legend(bbox_to_anchor = (1.5, 1))
    return f