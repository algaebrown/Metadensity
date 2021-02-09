import pyBigWig
from pybedtools import BedTool, create_interval_from_list
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import entropy
import pickle
import sys
import os
from collections import Counter, defaultdict


from .readdensity import ReadDensity
from .truncation import read_start_sites
from . import settings
transcript = BedTool(os.path.join(settings.root_dir, settings.transcript_fname))
gencode_feature = BedTool(os.path.join(settings.root_dir, settings.gencode_feature_fname))
print('Using: ', os.path.join(settings.root_dir, settings.transcript_fname))

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
        
class STAMP:
    """
    STAMP object to contain edit sites
    """
    def __init__(self):
        self.read_densities = {}
        self.name = ''
        self.uID = ''
        self.single_end = False
        self.rep_keys = []
        self.peaks = {} # edits for STAMP
        self.idr = None
        
class eCLIP:
    """ 
    eCLIP object to cotain .bigwig, .bam, .bed
    """
    def __init__(self, name='RBP', uID=None, single_end=True, rep_keys = ['rep1', 'rep2'], read_densities={}, peaks = {}, idr = None):
        # experimental variables
        self.name = name
        if uID:
            self.uID = uID
        else:
            self.uID = name
        self.rep_keys = rep_keys
        self.single_end = single_end

        # data
        self.read_densities = read_densities # bigwig, bam
        self.peaks = peaks # bed 
        self.idr = idr # bed
    
    @classmethod
    def from_series(cls, series, single_end = True):
        ''' build eCLIP object from series containing filenames for bam, bw, beds,
        series index must be :
        RBP for name; uid for uID.
        bam_0, bam_1, bam_control (or bam_control_0, bam_control_1);
        plus_0, plus_1, plus_control for pos.bw; minus_0, minus_1, minus_control for neg.bw;
        bed_0, bed_1 for individual peaks, idr for IDR peaks.
        _0, _1 is for replicate 1,2; you can extend to the # of reps you have.
        for control, if you have more than 1 SMInput, use bam_control_0, bam_control_1...etc.
        If only 1 SMInput, simply omit _0. use bam_control, plue_control...etc
        '''

        # autodetect SMInput structure
        if 'bam_control' in series.index:
            single_control = True
        else:
            single_control = False
        # autodetect replicate number
        n_reps = len([s for s in series.index if 'bam' in s and 'control' not in s]) # bam_1, bam_0
        assert n_reps > 0
        rep_keys = ['rep{}'.format(n+1) for n in range(n_reps)] # rep1 and rep2

        # read densities!
        read_densities = {}
        for i, rep in enumerate(rep_keys):
            pos = series['plus_{}'.format(i)]
            neg = series['minus_{}'.format(i)]
            bam = series['bam_{}'.format(i)]
            read_densities[rep] = ReadDensity(pos, neg, name = series['RBP'], bam = bam)
            if not single_control:
                pos = series['plus_control_{}'.format(i)]
                neg = series['minus_control_{}'.format(i)]
                bam = series['bam_control_{}'.format(i)]
                read_densities['ctrl{}'.format(i+1)] = ReadDensity(pos, neg, name = series['RBP'], bam = bam) # ctrl1, ctrl2 paired with rep1, rep2 ...
        if single_control:
            pos = series['plus_control']
            neg = series['minus_control']
            bam = series['bam_control']
            read_densities['ctrl'] = ReadDensity(pos, neg, name = series['RBP'], bam = bam)

        # peaks
        peaks = {}
        for i, rep in enumerate(rep_keys):
            peaks[rep] = BedTool(series['bed_{}'.format(i)])
        
        # idr
        idr = BedTool(series['idr'])

        return cls(name=series['RBP'], uID=series['uid'], single_end=single_end, rep_keys = rep_keys, read_densities=read_densities, peaks = peaks, idr = idr)

        

    
    def find_idr_transcript(self, genome_coord = transcript):
        """ find transcript containing at least 1 IDR peak
        Kwargs:
            genome_coord: (BedTool) transcript coordintae. default are canonical transcripts.
        """
        self.idr_transcript = genome_coord.intersect(self.idr, s=True, wa = True, u = True).saveas()

    ##################### FOR ENCODE DATA USE ONLY ###################################
        
    
        
        
### Metadensity

class Meta:
    ''' superclass of to inherit from '''
    def __init__(self, eCLIP, name, sample_no = 200, metagenes = None, transcripts = None):
        
        self.eCLIP = eCLIP
        
        self.name = name
        
        
        # TODO when initializing the exons are empty for some classes; but called using child behave normally.??
        if metagenes:
            self.metagene = metagenes # precomupted 
        elif transcripts:
            self.transcript = transcripts # BedTools object
            self.build_metagene(sample_no = sample_no)
        else:
            # by default use IDR transcripts
            self.eCLIP.find_idr_transcript()
            self.transcript = self.eCLIP.idr_transcript
            self.build_metagene(sample_no = sample_no)
        
        # get feature names
        self.featnames = list(self.metagene.values())[0].featnames
    
    def build_metagene(self, sample_no = 200):
        ''' Use function `Build_many_metagenes()` to get regions in UTR, intron, exon for each transcripts in self.idr_transcript and self.no_peak; store metagene in self.idr_metagene/self.neg_metagene'''
        tids = [s.attrs['transcript_id'] for s in self.transcript]
        self.metagene = Build_many_metagene(tids, sample_no = sample_no)
    def get_feature_density_array(self, feature, target_len, align, pad_value = np.nan, use_quantile = False, use_truncation = False):
        ''' align features by padding np.nan to the 5' or 3' end '''
        
        metagenes = self.metagene.values()
        for rep in self.eCLIP.rep_keys:
            rep_densities = []
            for m in metagenes:          
                        
                if use_quantile:
                    den = m.qdensities[self.eCLIP.uID][rep][feature]
                if use_truncation:
                    den = m.truncations[self.eCLIP.uID][rep][feature]
                else:
                    den = m.densities[self.eCLIP.uID][rep][feature]
            
            
                              
                # see if multi-feature average is called
                if type(den) == tuple:
                    if align == 'left':
                        rep_densities.append(
                            trim_or_pad(den[0], target_len, align, pad_value)
                        )
                        
                    elif align == 'right':
                        rep_densities.append(
                            trim_or_pad(den[1], target_len, align, pad_value)
                        )
                        
                else: # only 1 entry, no such align problem!
                    rep_densities.append(
                        trim_or_pad(den, target_len, align, pad_value)
                    )
                 
        
            # store differntly
            if use_quantile:
                self.qdensity_array[feature,align,rep]= np.stack(rep_densities)
            if use_truncation:
                self.truncate_array[feature,align,rep]= np.stack(rep_densities)
            else:
                self.density_array[feature,align,rep]= np.stack(rep_densities)
    def get_density_array(self, use_quantile = False, use_truncation = False):
        ''' extract metadensity from each metagene, zero pad or trim to get np array '''
        
        
        if use_quantile:
            self.qdensity_array = {}
        if use_truncation:
            self.truncate_array = {}
        else:
            self.density_array = {}
        
        for feature in self.featnames:
            for align in ['left', 'right']:
                flen = settings.feat_len[feature]
                self.get_feature_density_array(feature, flen , align, use_quantile = use_quantile, use_truncation = use_truncation)
    def concat_density_array(self, rep = 'rep1', quantile = False):
        ''' return concatenated density array by sequence of featnames '''
        # TODO make compatible to truncation 
        if quantile:
            concat = np.concatenate([self.qdensity_array[feat, align, rep] 
                        for feat in self.featnames for align in ['left', 'right']],axis = 1)
        else:
            
            concat = np.concatenate([self.density_array[feat, align, rep] 
                        for feat in self.featnames for align in ['left', 'right']],axis = 1)
        return concat
    def scale_density_array(self, method = 'norm_by_sum', quantile = False):
        
        ''' scale density array by only the features considered '''
        # TODO make compatible to truncation 
        self.scaled_density_array = {}
        if quantile:
            
            denarray = self.qdensity_array
        else:
            
            denarray = self.density_array
            
        for rep in self.eCLIP.rep_keys:
        
            rep_concat = self.concat_density_array(rep = rep, quantile = quantile)
            
        
            if method == 'norm_by_sum':
                scale = np.nansum(rep_concat, axis = 1)[:, None]
                
            
            if method == 'max_scalar':
                scale = np.nanmax(rep_concat, axis = 1)[:, None]
                
                      
            for align in ['left', 'right']:
                for feature in self.featnames:
                    self.scaled_density_array[feature,align,rep] = denarray[feature,align,rep]/scale

class Metatruncate(Meta):
    ''' Metagene for 5' read start '''
    def __init__(self, eCLIP, name, sample_no = 200, background_method = 'subtract', normalize = True, metagenes = None, transcripts = None, smooth = False):
        # generate metagene coords
        super().__init__(eCLIP, name, sample_no, metagenes = metagenes, transcripts = transcripts)

        # automatically run
        self.get_truncation(background_method = background_method, normalize = normalize, smooth = smooth)
        self.background_method = background_method
        self.normalize = normalize
    
    
    def get_truncation(self, background_method = 'subtract', normalize = True, smooth = False):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object
        background_method: 'subtract'; 'subtract normal'; 'relative information' How to deal with IP and Input. None for doing nothing.
        normalize: True to use % of edits in position; False to use raw background subtract signals
        '''
        _ = [m.get_average_feature(self.eCLIP, background_method = background_method, normalize = normalize, truncate = True, smooth = smooth) for m in self.metagene.values()]
    
    
    
class Metadensity(Meta):
    ''' Metadensity can be created from eCLIP and a set of transcripts'''
    def __init__(self, eCLIP, name, sample_no = 200, background_method = 'subtract', normalize = True, metagenes = None, transcripts = None, smooth = False):
        # generate metagene coords
        super().__init__(eCLIP, name, sample_no, metagenes = metagenes, transcripts = transcripts)

        # automatically run
        self.get_metadensity(background_method = background_method, normalize = normalize, smooth = smooth)
        self.background_method = background_method
        self.normalize = normalize
    
        
    def get_metadensity(self, background_method = 'subtract', normalize = True, smooth = False):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object'''
        _ = [m.get_average_feature(self.eCLIP, background_method = background_method, normalize = normalize, truncate = False, smooth = smooth) for m in self.metagene.values()]
    
    
    
    
    

        
        

########################## Metagene class #####################################

class Metagene:
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        self.id = esnt # can be other depending on stuffs
        self.type = transcript_type
        self.chrom = chro
        self.start = start
        self.stop = end
        self.strand = strand

        if not features:
            self.features = defaultdict(set)
        else:
            self.features = features
        
        self.featnames = ['first_exon', 'exon', 'intron', 'branchpoint', 'branchpoint_pred','last_exon']
        
        

        # TODO record raw value before pooling, save I/O time
        self.coverage = {} # key: eCLIP uID, values: coverage RAW from 5' end to 3' end
        self.sites = {} # key: eCLIP uID, values: RAW truncation count
        
        # pooled metagene density from such
        self.densities = {} # key: eCLIP uID, values: metadensity
        self.qdensities  = {}
        self.truncations = {}
    
    @classmethod
    def from_dict(cls, coord_dict):
        return cls(coord_dict['id'],
        coord_dict['chrom'], coord_dict['start'], coord_dict['end'], coord_dict['strand'],
         coord_dict['type'], features = coord_dict['features'])
    ################################## about annotation ###########################################
    
    
    def first_last_exon(self):
        ''' find first, last exon and cds '''
        # seperate first and last exon
        if len(self.features['exon']) != 0:
            min_start = min([e[0] for e in list(self.features['exon'])])
            max_start = max([e[1] for e in list(self.features['exon'])])
        
            if self.strand == '+':
                self.features['first_exon'] = set([e for e in list(self.features['exon']) if e[0] == min_start])
                self.features['last_exon'] = set([e for e in list(self.features['exon']) if e[1] == max_start])
            else:
                self.features['last_exon'] = set([e for e in list(self.features['exon']) if e[0] == min_start])
                self.features['first_exon'] = set([e for e in list(self.features['exon']) if e[1] == max_start])
        
            self.features['exon'] = self.features['exon'] - self.features['first_exon'] - self.features['last_exon']
        else:
            self.features['first_exon'] = set()
            self.features['last_exon'] = set()
        if len(self.features['CDS']) != 0:
            min_start = min([e[0] for e in list(self.features['CDS'])])
            max_start = max([e[1] for e in list(self.features['CDS'])])
        
            if self.strand == '+':
                self.features['first_CDS'] = set([e for e in list(self.features['CDS']) if e[0] == min_start])
                self.features['last_CDS'] = set([e for e in list(self.features['CDS']) if e[1] == max_start])
            else:
                self.features['last_CDS'] = set([e for e in list(self.features['CDS']) if e[0] == min_start])
                self.features['first_CDS'] = set([e for e in list(self.features['CDS']) if e[1] == max_start])
        
            self.features['CDS'] = self.features['CDS'] - self.features['first_CDS'] - self.features['last_CDS']
        else:
            self.features['first_CDS'] = set()
            m=self.features['last_CDS'] = set()

        
    def five_or_three_utr(self):
        ''' to distinguish five or three utr '''
        
        if len(self.features['UTR']) != 0:
            
            # seperate first and last utr
            min_start = min([e[0] for e in list(self.features['UTR'])])
            max_start = max([e[1] for e in list(self.features['UTR'])])

            if self.strand == '+':
                self.features['five_utr'] = set([e for e in list(self.features['UTR']) if e[0] == min_start])
                self.features['three_prime_UTR'] = set([e for e in list(self.features['UTR']) if e[1] == max_start])
            else:
            
                self.features['five_utr'] = set([e for e in list(self.features['UTR']) if e[1] == max_start])
                self.features['three_prime_UTR'] = set([e for e in list(self.features['UTR']) if e[0] == min_start])
    
    ################################## about raw values ###########################################
    def truncation_count(self, eCLIP, rep, feature):
        ''' return a list of read start count for each position in feature
        ex: sites return [10,11,12,13]; feature = (10,15); return [1,1,1,1,0] count from 5' to 3'
        '''
        if isinstance(eCLIP, STAMP):
            s = '{}\t{}\t{}\t{}\t{}\t{}'.format(self.chrom, str(feature[0]), str(feature[1]), '.', '.', self.strand)
            bed = BedTool(s, from_string=True).saveas()
            sites_bed = eCLIP.peaks[rep].intersect(bed, s = True).saveas()
            sites = [int(s.start) for s in sites_bed]
            
        else:
            
            sites = read_start_sites(eCLIP.read_densities[rep].bam, chrom=self.chrom, start=feature[0], end=feature[1], strand=self.strand, single_end = eCLIP.single_end)
        
        site_count = Counter(sites)
        
        pos_count = []
        if self.strand == '+':
            for pos in range(feature[0], feature[1]+1):
                pos_count.append(site_count[pos])
        if self.strand == '-':
            for pos in range(feature[1], feature[0]-1, -1): # 15, 14, 13, 12, 11, 10
                pos_count.append(site_count[pos])
        return np.array(pos_count)  
    
    def coverage_value(self, eCLIP, rep, feature):
        ''' return raw RPM value from bigwig files for each position in feature 
        from 5' to 3' end '''
        rpm = np.nan_to_num(eCLIP.read_densities[rep].values(self.chrom, self.start, self.stop, self.strand),0)
        if self.strand == '-':
            rpm = -rpm
        
        return rpm
    
    def get_raw_value(self, eCLIP, truncate = False):
        ''' fetch raw coverage/truncation value for whole gene 
        write in self.coverage or self.sites'''
        
        if truncate:
            # check if it has been computed (save I/O time)
            if eCLIP.uID in self.sites.keys(): # already computed
                return 1 # terminate fetching
            else:
            
                self.sites[eCLIP.uID] = {}
                save = self.sites
        else:
            if eCLIP.uID in self.coverage.keys():
                return 1
            else:
                self.coverage[eCLIP.uID] = {}
                save = self.coverage
        
        # for ips and controls
        for rep in eCLIP.read_densities.keys():
            
            if truncate:
                save[eCLIP.uID][rep] = self.truncation_count(eCLIP, rep, (self.start, self.stop)) # save list from 5' to 3'
            else:
                save[eCLIP.uID][rep] = self.coverage_value(eCLIP, rep, (self.start, self.stop))
    ################################## about background removal ###########################################
    def remove_background(self, eCLIP, rep, method = 'subtract', truncate = False):
        ''' remove background by method
        method = 'subtract': subtract IP from Input
        method = 'subtract normal': subtract IP distribution from Input distribution
        method = 'relative information': take relative entropy
        method = 'get null': get input density
        '''
        if truncate:
            save = self.sites
        else:
            save = self.coverage
        
        # determine how many control
        if 'ctrl' in save[eCLIP.uID].keys():
            control_raw = save[eCLIP.uID]['ctrl']
        else:
            control_raw = save[eCLIP.uID]['ctrl'+rep.replace('rep','')] # ctrl1
        
        raw = save[eCLIP.uID][rep]
        
        # perform background removal
        if method == 'subtract':
            values = raw - control_raw
            # replace negative value with 0
            values[values < 0] = 0
        if method == 'subtract normal':
            if np.sum(raw) == 0:
                values = np.zeros(raw.shape)
            elif np.sum(control_raw) == 0:
                values = raw/np.sum(raw)
            else:
                values = raw/np.sum(raw) - control_raw/np.sum(control_raw)
            # replace negative value with 0
            values[values < 0] = 0
        if method == 'relative information':
            # handle 0, add pseudocount
            
            if truncate:
                ip_pseudocount = 1 #TODO this is not compatible with density
                control_pseudocount = 1
            else:
                non_zero = raw[np.where(raw>0)]
                if len(non_zero) == 0:
                    ip_pseudocount = 1 # doesn't matter, uniform distribution
                else:
                    ip_pseudocount = np.min(non_zero)

                non_zero = control_raw[np.where(control_raw>0)]
                if len(non_zero) == 0:
                    control_pseudocount = 1 # doesn't matter, uniform distribution
                else:
                    control_pseudocount = np.min(non_zero)
            
            raw_added = raw + ip_pseudocount
            control_added = control_raw + control_pseudocount

            # get fraction, pi and qi
            ip_norm = raw_added/np.sum(raw_added) # pi
            input_norm = control_added/np.sum(control_added) # qi
            
            


            values = np.multiply(
                ip_norm,
                np.log2(np.divide(ip_norm, input_norm))
                ) # pi * log2(pi/qi)
        if method == None:
            values = raw # do nothing with background
        
        return values   
    
    ################################## about axis ###########################################    
    def to_relative_axis(self, feature):
        ''' self.start = 1000, self.stop = 5000, self.strand = -, feature = (1000,1005)
        return (5000-1005, 5000-1000) = (3995, 4000)
        if strand +, (1000-1000, 1005-1000) '''
        if self.strand == '-':
            relative_feature = (self.stop - feature[1], self.stop - feature[0])
        else:
            relative_feature = (feature[0]-self.start, feature[1] - self.start)
        return relative_feature
    
    ################################## about annotation ###########################################
    def multi_feature_avg(self, feature, values, align = 'left', max_len = None):
        '''
        average and align (zero padding) multiple intron/exon
        feature: list of intervals for multiple exons/cds/intron; ex: [(1,100), (2,950)]
        values: 5' to 3' processed value (remove background, normalize) of the whole gene
        align: align feature of different length to the left or right
        return: average feature (np.array)
        '''
        all_feature_values = []
        # feature length
        if max_len == None:
            if type(list(feature)[0]) == int:
                max_len = 20
            else:
                max_len = max([f[1]-f[0] for f in feature])
        
        for f in feature:
            if type(f) == int:
                f = (f-10, f+10) # single point feature

            axis = self.to_relative_axis(f)
            feature_values = values[axis[0]:axis[1]] # window around
            all_feature_values.append(trim_or_pad(feature_values, max_len, pad_value = np.nan, align = align)) # pad to the same length then append
        
        
        
        feature_average  = np.nanmean(np.stack(all_feature_values), axis = 0)
        return feature_average
        
    

        
    def concat_density(self, uID, rep, truncation = False):
        ''' since some feature would have left/right, we need a special function to return concat density '''
        if truncation:
            left_densities = [feat[0] if type(feat) == tuple else feat for feat in self.truncations[uID][rep].values()]
        else:
            
            left_densities = [feat[0] if type(feat) == tuple else feat for feat in self.densities[uID][rep].values()]
        return np.concatenate(left_densities)
        

    
    def get_average_feature(self, eCLIP, background_method = 'subtract', normalize = True, truncate = True, smooth = False):
        ''' fetching values for each feature, also averaging them if needed 
        write values in self.truncate or self.densities, depending on truncate'''
        # fetch raw sites for all reps and ctrl
        self.get_raw_value(eCLIP, truncate = truncate)

        # save to different attributes
        if truncate:
            to_save = self.truncations
            raw = self.sites
        else:
            to_save = self.densities
            raw = self.coverage
        
        

        # initialize empty dictionary to store result
        to_save[eCLIP.uID] = {}
        def nan_array():
            minus1 = np.empty(1)    
            minus1[:] = np.nan
            return minus1
        for rep in eCLIP.rep_keys:
            to_save[eCLIP.uID][rep] = defaultdict(nan_array)

            # deal with background
            if background_method == 'get null':
                # just to get input distribution
                if 'ctrl' in raw[eCLIP.uID].keys():
                    values = raw[eCLIP.uID]['ctrl']
                else:
                    values = raw[eCLIP.uID]['ctrl{}'.format(rep.replace('rep', ''))]
            elif background_method: # if background method is not None, which means not dealing with background
                values = self.remove_background(eCLIP, rep, method = background_method, truncate = truncate)
            else:
                values = raw[eCLIP.uID][rep]
            # smooth
            if smooth:
                values = gaussian_smooth(values)

            # deal with normalization
            if normalize:
                values = values/np.sum(values)
                # TODO: add quantile for density
            else:
                pass
            
            # getting features:
            for feature, fname in zip(self.features.values(), self.features.keys()):
            
                if len(feature) == 0: # no such feature
                    # np.empty does not always yield nan
                    
                    to_save[eCLIP.uID][rep][fname] = nan_array()
                elif len(feature) == 1:
                    feature = list(feature)[0]
                    if type(feature) == int: # single point feature
                        feature = (feature-10, feature + 10) # extend a window for it
                    axis = self.to_relative_axis(feature)
                
                    
                    minus1 = values[axis[0]: axis[1]]
                    to_save[eCLIP.uID][rep][fname] = minus1
                else:
                    if type(feature) == int: # single point feature
                        # no issue with alignment, single point feature always have the same length
                        left_mean = self.multi_feature_avg(feature, values, align = 'left') 
                        to_save[eCLIP.uID][rep][fname] = (left_mean, right_mean)
                    else:
                        left_mean = self.multi_feature_avg(feature, values, align = 'left') 
                        right_mean = self.multi_feature_avg(feature, values, align = 'right')
                

                        to_save[eCLIP.uID][rep][fname] = (left_mean, right_mean)
        
 
    def quantile_metadensity(self, q = 40):
        
        ''' DEPRECATED quantile metadensity to equal sized bins at transcript level'''
        
        for uID in self.densities.keys():
            
            self.qdensities[uID] = {}
            for rep in self.densities[uID].keys():
                self.qdensities[uID][rep] = {}
                
                
                concat = self.concat_density(uID, rep)
                try:
                    # quantile at a transcript level
                    quantile, bins = pd.qcut(concat, q, retbins = True,duplicates = 'drop')
                    for feat in self.featnames:
                        density = self.densities[uID][rep][feat]
                        if type(density) == tuple:
                            # accomodate mutliple intron/exon; has left right property
                            self.qdensities[uID][rep][feat] = tuple([np.digitize(d,bins) for d in density])
                            
                        else:
                            self.qdensities[uID][rep][feat]  = np.digitize(density,bins)
                except:
                    
                    print(self.ensembl_id, 'unique values: ', len(np.unique(concat)))
                    for feat in self.featnames:
                        self.qdensities[uID][rep][feat]  = np.empty(1)
                        self.qdensities[uID][rep][feat][:] = np.nan # np.empty can generate some weirdly large number
    
    def quantile_metadensity_raw(self, eCLIP, q = 40):
        ''' DEPRECATED quantilize metadensity uning raw density (not pooled, not averaged)'''
        uID = eCLIP.uID
        
        # initialize empty dictionary
        self.qdensities[uID] = {}
        feature_den1 = []
        feature_den2 = []
        
        # get bins
        result_dict= remove_background(eCLIP, self.chrom, self.start, self.stop, self.strand)
        transcript1 = result_dict['rep1']
        transcript2 = result_dict['rep2']
        quantile1, bins1 = pd.qcut(transcript1, q, retbins = True,duplicates = 'drop')
        quantile2, bins2 = pd.qcut(transcript2, q, retbins = True,duplicates = 'drop')
        
        # warn if don't reach q bins
        if len(bins1) < q+1 or len(bins2) < q+1:
            import warnings
            warnings.warn('Fail to reach {} bins, transcript {} has only {} {} unique values, {} {} non-nan values'.format
                         (q, self.ensembl_id,
                          len(np.unique(transcript1)), len(np.unique(transcript2)),
                          len(transcript1) - np.sum(np.isnan(transcript1)), len(transcript2) - np.sum(np.isnan(transcript2))
                         )
                         )
        
        # get average density per feature
        for feature in self.features.values():
            if len(feature) == 0: # no such feature
                minus1 = np.empty(1) # np.empty does not always yield nan
                minus2 = np.empty(1)
                minus1[:] = np.nan
                minus2[:] = np.nan
                
            elif len(feature) == 1:
                feature = list(feature)[0]
                result_dict = remove_background(eCLIP, self.chrom, feature[0], feature[1], self.strand)
                minus1 = result_dict['rep1']
                minus2 = result_dict['rep2']
                minus1 = np.array(np.digitize(minus1, bins1)) #  np.digitize(self.densities[uID][rep][feat],bins)
                minus2 = np.array(np.digitize(minus2, bins2))
            else:
                minus1_l, minus2_l = self.multi_feature_avg(feature, eCLIP, quantile = True, bins1 = bins1, bins2 = bins2, align = 'left') 
                minus1_r, minus2_r = self.multi_feature_avg(feature, eCLIP, quantile = True, bins1 = bins1, bins2 = bins2, align = 'right') 
                
                minus1 = (minus1_l, minus1_r)
                minus2 = (minus2_l, minus2_r)
            feature_den1.append(minus1)
            feature_den2.append(minus2)
        
        # No need to normalize after quantilization
        
        n1 = feature_den1
        n2 = feature_den2
        
    
    
        fnames = self.features.keys()
        
        self.qdensities[eCLIP.uID] = {
                                    'rep1':dict(zip(fnames, n1)),
                                    'rep2':dict(zip(fnames, n2))
                                    }
class Protein_coding_gene(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
        self.featnames+= ['five_prime_UTR','first_CDS', 'CDS', 'last_CDS', 'three_prime_UTR']
        

class vaultRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class lncRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class RNU(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class circRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class YRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class snoRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
        self.featnames += ['5p_duplex', '3p_duplex', 'mature'] # TODO the boxes
        
        
class rRNA(Metagene):
    # TODO get from repetitive pipeline
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class tRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
        self.featnames += ['D_loop', 'T_loop', 'Anticodon_loop', 'Variable_loop']

class miRNA(Metagene):
    ''' a class for miRNA's features specifically '''
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
        self.featnames += ['5p_duplex', '3p_duplex', 'mature']

    def find_hairpin(self):
        if len(self.fivep_duplex)>0 and len(self.threep_duplex)>0:
            if self.strand == '+':
                # 5' will have smaller number
                start = list(self.fivep_duplex)[0][1] # end of 5p duplex
                end = list(self.threep_duplex)[0][0] # start of 3p duplex
            if self.strand == '-':
                start = list(self.threep_duplex)[0][1] # end of 5p duplex
                end = list(self.fivep_duplex)[0][0] # start of 3p duplex
            
            self.hairpin.update([(start, end)])
            # TODO download structure from miRbase to find real hairpins



    
def Build_many_metagene(id_list, sample_no = 200, transcript_type = 'protein_coding'):
    ''' Create Metagene object for regions for key transcript 
    hg19 does not have 'three_prime_UTR and five_prime_UTR which is annoying'''
    # Build empty metagene object
    all_metagene = {}
   
    if len(id_list)< sample_no:
        ids = id_list
        
    else: 
        ids = id_list[:sample_no]

    # load pre-parsed gencode coords and branchpoint
    
    if transcript_type == 'miRNA':
        coord_dicts = pickle.load(open(settings.mir, 'rb'))
    elif transcript_type == 'snoRNA':
        coord_dicts = pickle.load(open(settings.sno, 'rb'))
    else:
        print('Using:',settings.gencode)
        coord_dicts = pickle.load(open(settings.gencode, 'rb'))
    
    for i in ids:
        d = coord_dicts[i]
        if d['type'] == transcript_type:
            d['features'] = defaultdict(set, d['features']) # convert to default dict in case feature doesn't exist
            if transcript_type == 'protein_coding':
                all_metagene[i] = Protein_coding_gene.from_dict(d)# TODO: miR, snoR will need another dictionary
            elif transcript_type == 'miRNA':
                all_metagene[i] = miRNA.from_dict(d)
            else:
                all_metagene[i] = Metagene.from_dict(d)

    print('Done building metagene')
    
    return all_metagene

    




########################################## calculate meta-density ##########################################

def trim_or_pad(density, target_length, align = 'left', pad_value = 0):
    ''' make density your target length by trimming or padding'''
    if len(density) == target_length:
        return density
    elif len(density) > target_length:
        if align == 'left':
            return density[:target_length]
        if align == 'right':
            return density[-target_length:]
    else:
        discrepency = target_length - len(density)
        if align == 'left':
            density = np.array(list(density) + [pad_value]*discrepency)
            
            
        if align == 'right':
            density = np.array([pad_value]*discrepency +list(density))
        return density


######################################## Probability Distribution #########################################
def bind_strength_discretize(density_array, bins = 8, ymax = 0.03):
    ''' given density array(sample * length of feature), count discretization'''
    pseudocount = [np.histogram(d[~np.isnan(d)], range = (0, ymax), bins = bins)[0] + 1 for d in density_array.T] # pseudocount 
    den_prob = np.stack([p/np.sum(p) for p in pseudocount]).T
    
    return den_prob

def pos_spec_bind_strength(eCLIP, peak_max = 0.03, bins = 20, use_quantile = False, use_scaled = False):
    ''' return probability distribution on \"binding stregth/normalized peak height\" for each position'''
    bind_strength = {}
    
    if use_quantile:
        den_array = eCLIP.qdensity_array
    if use_scaled:
        den_array = eCLIP.scaled_density_array
    
    else:
        den_array = eCLIP.density_array
    for k in den_array.keys():
        bind_strength[k] = bind_strength_discretize(den_array[k], bins = bins, ymax = peak_max)
    return bind_strength
######################################### Relative entropy ##########################################


def density_array_entropy(f, fn):
    '''
    relative entropy for each position, f = positive binding density; fn = negative binding examples's density
    '''
    return np.array([entropy(pk = f[1:, pos], qk = fn[1:, pos]) for pos in range(f.shape[1])])

