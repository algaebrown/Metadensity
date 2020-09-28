import pyBigWig
import pandas as pd
import sys
sys.path.append('/home/hsher/rbp-maps/maps/')

from density.ReadDensity import ReadDensity
from truncation import read_start_sites
import os
from pybedtools import BedTool, create_interval_from_list
import matplotlib.pyplot as plt
import numpy as np
basedir = '/home/hsher/seqdata/eclip_raw/'
from scipy.stats import entropy
from collections import Counter
from config import settings
transcript = BedTool(os.path.join(settings.root_dir, settings.transcript_fname))
gencode_feature = BedTool(os.path.join(settings.root_dir, settings.gencode_feature_fname))

featnames = ['five_utr', 'first_exon', 'exon', 'intron', 'last_exon','three_utr', 'CDS', 'first_CDS', 'last_CDS']

def make_density(series, basedir = basedir):
    ''' Generate 3 ReadDensity Object from pd.Series from encode_data_id.pickle'''
    all_den = []
    for types in ['0','1','control']: # rep1, rep2, control
        neg = basedir + series['minus_'+types].values[0]
        pos = basedir + series['plus_'+types].values[0]
        bam = basedir + series['bam_'+types].values[0]
        
        density = ReadDensity(pos, neg, bam = bam, name = str(series['RBP']))
        all_den.append(density)
    return all_den[0], all_den[1], all_den[2]
        
class STAMP:
    def __init__(self):
        self.read_densities = {}
        self.name = ''
        self.uID = ''
        self.rep_keys = []
        self.peaks = {} # edits for STAMP
        self.idr = None
        self.single_end = False
class eCLIP:
    def __init__(self):
        self.read_densities = {}
        self.name = ''
        self.uID = ''
        self.rep_keys = []
        self.peaks = {}
        self.idr = None
        self.single_end = False
        
    def build_ENCODE(self, series):
        ''' from the above ENCODE data dataframe'''
        rep1, rep2, ctrl = make_density(series)
        self.read_densities['rep1'] = rep1
        self.read_densities['rep2'] = rep2
        self.read_densities['ctrl'] = ctrl
        self.uID = series['uID'].values[0]
        self.name = series['RBP'].values[0]

        self.rep_keys = ['rep1', 'rep2']

        # automatically add peaks
        self.add_peaks()
    def build_YTHDF(self):
        ''' 3 control, 3 IP structure '''
        # to do here
        pass

    def add_peaks(self,
                  idr_path = '/home/hsher/seqdata/eclip_bed/sorted/',
                  indv_path = '/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20180205/hg38/'):
        ''' add peaks to eCLIP object'''
        idr_suffix = '.01v02.IDR.out.0102merged.bed.blacklist_removed.bed.narrowPeak.bed'
        self.idr = BedTool(idr_path + self.uID + idr_suffix)
        
        indv_suffix = '{0}_0{1}.basedon_{0}_0{1}.peaks.l2inputnormnew.bed.compressed.bed.blacklist_removed.bed'
        self.peak1 = BedTool(indv_path + indv_suffix.format(self.uID, 1))
        self.peak2 = BedTool(indv_path + indv_suffix.format(self.uID, 2))
    def find_idr_transcript(self, genome_coord = transcript):
        ''' find positive example for generateing data'''
        self.idr_transcript = genome_coord.intersect(self.idr, s=True, wa = True, u = True).saveas()
    def find_negative_example(self, genome_coord = transcript):
        ''' find negative regions as no IDR neither individual peaks'''
        self.no_peak = genome_coord.intersect(self.idr, s=True, v=True).intersect(self.peak1, s = True, v =True).intersect(self.peak2, s = True, v =True).saveas()
    
    
        
    def RBP_centric_approach(self, series, sample_no = 200):
        ''' create eCLIP object, get peaks, get examples, metagene and metadensity with one function. return eCLIP object'''
        self.build_ENCODE(series)
        
        print('finding negative/positive examples')
        self.find_idr_transcript()
        self.find_negative_example()
        
### Metadensity

class Meta:
    ''' superclass of to inherit from '''
    def __init__(self, eCLIP, transcripts, name, sample_no = 200, metagenes = None):
        self.eCLIP = eCLIP
        self.transcript = transcripts # BedTools object
        self.name = name
        
        # TODO when initializing the exons are empty for some classes; but called using child behave normally.??
        if metagenes:
            self.metagene = metagenes # precomupted 
        else:
            self.build_metagene(sample_no = sample_no)
    
    def build_metagene(self, sample_no = 200):
        ''' Use function `Build_many_metagenes()` to get regions in UTR, intron, exon for each transcripts in self.idr_transcript and self.no_peak; store metagene in self.idr_metagene/self.neg_metagene'''
        self.metagene = Build_many_metagene(self.transcript, sample_no = sample_no)
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
    def get_density_array(self, five_utr_len=100, three_utr_len=150, intron_len = 1500, exon_len = 150, cds_len = 150, use_quantile = False, use_truncation = False):
        ''' extract metadensity from each metagene, zero pad or trim to get np array '''
        
        
        if use_quantile:
            self.qdensity_array = {}
        if use_truncation:
            self.truncate_array = {}
        else:
            self.density_array = {}
        
        for feature, l in zip(featnames, [five_utr_len, exon_len,exon_len, intron_len, exon_len,three_utr_len, cds_len, cds_len, cds_len]):
            for align in ['left', 'right']:
                
                self.get_feature_density_array(feature, l, align, use_quantile = use_quantile, use_truncation = use_truncation)
    def concat_density_array(self, rep = 'rep1', quantile = False):
        ''' return concatenated density array by sequence of featnames '''
        # TODO make compatible to truncation 
        if quantile:
            concat = np.concatenate([self.qdensity_array[feat, align, rep] 
                        for feat in featnames for align in ['left', 'right']],axis = 1)
        else:
            
            concat = np.concatenate([self.density_array[feat, align, rep] 
                        for feat in featnames for align in ['left', 'right']],axis = 1)
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
                for feature in featnames:
                    self.scaled_density_array[feature,align,rep] = denarray[feature,align,rep]/scale

class Metatruncate(Meta):
    ''' Metagene for 5' read start '''
    def __init__(self, eCLIP, transcripts, name, sample_no = 200, background_method = 'subtract', normalize = True, metagenes = None):
        # generate metagene coords
        super().__init__(eCLIP, transcripts, name, sample_no, metagenes = metagenes)

        # automatically run
        self.get_truncation(background_method = background_method, normalize = normalize)
        self.background_method = background_method
        self.normalize = normalize
    
    
    def get_truncation(self, background_method = 'subtract', normalize = True):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object
        background_method: 'subtract'; 'subtract normal'; 'relative information' How to deal with IP and Input. None for doing nothing.
        normalize: True to use % of edits in position; False to use raw background subtract signals
        '''
        _ = [m.get_average_feature(self.eCLIP, background_method = background_method, normalize = normalize, truncate = True) for m in self.metagene.values()]
    
    
    
class Metadensity(Meta):
    ''' Metadensity can be created from eCLIP and a set of transcripts'''
    def __init__(self, eCLIP, transcripts, name, sample_no = 200, background_method = 'subtract', normalize = True, metagenes = None):
        # generate metagene coords
        super().__init__(eCLIP, transcripts, name, sample_no, metagenes = metagenes)

        # automatically run
        self.get_metadensity(background_method = background_method, normalize = normalize)
        self.background_method = background_method
        self.normalize = normalize
    
        
    def get_metadensity(self, background_method = 'subtract', normalize = True):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object'''
        _ = [m.get_average_feature(self.eCLIP, background_method = background_method, normalize = normalize, truncate = False) for m in self.metagene.values()]
    
    
    
    
    

        
        

########################## Metagene class #####################################

class Metagene:
    def __init__(self, esnt, chro, start, end, strand, transcript_type):
        self.ensembl_id = esnt
        self.type = transcript_type
        self.chrom = chro
        self.start = start
        self.stop = end
        self.strand = strand
        self.five_utr = set()
        self.three_utr = set()
        self.intron = set()
        self.exon = set()
        self.utr = set() # to accomodate hg19
        self.cds = set() 

        # TODO record raw value before pooling, save I/O time
        self.coverage = {} # key: eCLIP uID, values: coverage RAW from 5' end to 3' end
        self.sites = {} # key: eCLIP uID, values: RAW truncation count
        
        # pooled metagene density from such
        self.densities = {} # key: eCLIP uID, values: metadensity
        self.qdensities  = {}
        self.truncations = {}
    
    ################################## about annotation ###########################################
    def first_last_exon(self):
        ''' find first, last exon and cds '''
        # seperate first and last exon
        if len(self.exon) != 0:
            min_start = min([e[0] for e in list(self.exon)])
            max_start = max([e[1] for e in list(self.exon)])
        
            if self.strand == '+':
                self.first_exon = set([e for e in list(self.exon) if e[0] == min_start])
                self.last_exon = set([e for e in list(self.exon) if e[1] == max_start])
            else:
                self.last_exon = set([e for e in list(self.exon) if e[0] == min_start])
                self.first_exon = set([e for e in list(self.exon) if e[1] == max_start])
        
            self.exon = self.exon - self.first_exon - self.last_exon
        else:
            self.first_exon = set()
            self.last_exon = set()
        if len(self.cds) != 0:
            min_start = min([e[0] for e in list(self.cds)])
            max_start = max([e[1] for e in list(self.cds)])
        
            if self.strand == '+':
                self.first_cds = set([e for e in list(self.cds) if e[0] == min_start])
                self.last_cds = set([e for e in list(self.cds) if e[1] == max_start])
            else:
                self.last_cds = set([e for e in list(self.cds) if e[0] == min_start])
                self.first_cds = set([e for e in list(self.cds) if e[1] == max_start])
        
            self.cds = self.cds - self.first_cds - self.last_cds
        else:
            self.first_cds = set()
            m=self.last_cds = set()
    def five_or_three_utr(self):
        ''' to distinguish five or three utr '''
        
        if len(self.utr) != 0:
            
            # seperate first and last utr
            min_start = min([e[0] for e in list(self.utr)])
            max_start = max([e[1] for e in list(self.utr)])

            if self.strand == '+':
                self.five_utr = set([e for e in list(self.utr) if e[0] == min_start])
                self.three_utr = set([e for e in list(self.utr) if e[1] == max_start])
            else:
            
                self.five_utr = set([e for e in list(self.utr) if e[1] == max_start])
                self.three_utr = set([e for e in list(self.utr) if e[0] == min_start])
    
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
        rpm = np.nan_to_num(eCLIP.read_densities[rep].values(chrom, start, stop, strand),0)
        if self.strand == '-':
            rpm = -rpm
        
        return rpm
    
    def get_raw_value(self, eCLIP, truncate = False):
        ''' fetch raw coverage/truncation value for whole gene 
        write in self.coverage or self.sites'''
        
        if truncate:
            self.sites[eCLIP.uID] = {}
            save = self.sites
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
                ip_pseudocount = np.min(non_zero) # assume the minimal non-zero density equals one read

                non_zero = raw[np.where(control_raw>0)]
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
            max_len = max([f[1]-f[0] for f in feature])
        
        for f in feature:
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
        
    def protein_coding_gene(self, eCLIP, method = 'subtract', truncate = False):
        ''' return 5' UTR -first_cds - CDS concat - last_cds - 3' UTR metagene'''
        pass
        # TODO
    def generic_RNA(self, eCLIP, method = 'subtract', truncate = False):
        ''' return First exon - intron - exon -last exon '''
        pass
        # TODO
    
    def get_average_feature(self, eCLIP, background_method = 'subtract', normalize = True, truncate = True):
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
        for rep in eCLIP.rep_keys:
            to_save[eCLIP.uID][rep] = {}

            # deal with background
            if background_method: # if background method is not None, which means not dealing with background
                values = self.remove_background(eCLIP, rep, method = background_method, truncate = truncate)
            else:
                values = raw[eCLIP.uID][rep]

            # deal with normalization
            if normalize:
                values = values/np.sum(values)
                # TODO: add quantile for density
            else:
                pass
            
            # getting features:
            for feature, fname in zip([self.five_utr, self.first_exon, self.exon, self.intron, self.last_exon, self.three_utr, self.cds, self.first_cds, self.last_cds], featnames):
            
                if len(feature) == 0: # no such feature
                    minus1 = np.empty(1) # np.empty does not always yield nan
                    minus1[:] = np.nan
                    to_save[eCLIP.uID][rep][fname] = minus1
                elif len(feature) == 1:
                    feature = list(feature)[0]
                    axis = self.to_relative_axis(feature)
                
                    
                    minus1 = values[axis[0]: axis[1]]
                    to_save[eCLIP.uID][rep][fname] = minus1
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
                    for feat in featnames:
                        density = self.densities[uID][rep][feat]
                        if type(density) == tuple:
                            # accomodate mutliple intron/exon; has left right property
                            self.qdensities[uID][rep][feat] = tuple([np.digitize(d,bins) for d in density])
                            
                        else:
                            self.qdensities[uID][rep][feat]  = np.digitize(density,bins)
                except:
                    
                    print(self.ensembl_id, 'unique values: ', len(np.unique(concat)))
                    for feat in featnames:
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
        for feature in [self.five_utr, self.first_exon, self.exon, self.intron, self.last_exon, self.three_utr]:
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
        
    
    
        fnames = featnames
        
        self.qdensities[eCLIP.uID] = {
                                    'rep1':dict(zip(fnames, n1)),
                                    'rep2':dict(zip(fnames, n2))
                                    }
        

    


def Build_many_metagene(key_transcript, gencode_feature = gencode_feature, sample_no = 200):
    ''' Create Metagene object for regions for key transcript 
    hg19 does not have 'three_prime_UTR and five_prime_UTR which is annoying'''
    # Build empty metagene object
    all_metagene = {}
   
    if len(key_transcript)< sample_no:
        samples = key_transcript
        
    else: 
        samples = key_transcript[:sample_no]
    
    
    
    for s in samples:
        enst = s.attrs['transcript_id']
        transcript_type = s.attrs['transcript_type']
        metagene = Metagene(enst, s.chrom, s.start, s.stop, s.strand, transcript_type)
        all_metagene[enst] = metagene
    
            
    trigger_hg19=False
    # write where is intron, exon...
    for i in gencode_feature.filter(lambda x: x.attrs['transcript_id'] in all_metagene.keys()):
        enst = i.attrs['transcript_id']
        feature_type = i.fields[2]
        
        if feature_type == 'transcript':
            all_metagene[enst].intron.update([(i.start, i.stop)])
        if feature_type == 'exon':
            all_metagene[enst].exon.update([(i.start, i.stop)])
        if feature_type == 'five_prime_UTR':
            all_metagene[enst].five_utr.update([(i.start, i.stop)])
        if feature_type == 'three_prime_UTR':
            all_metagene[enst].three_utr.update([(i.start, i.stop)])
        if feature_type == 'CDS': 
            all_metagene[enst].cds.update([(i.start, i.stop)])
        if feature_type == 'UTR': # hate it for hg 19
            all_metagene[enst].utr.update([(i.start, i.stop)])
            trigger_hg19=True
    if trigger_hg19:
        # need to distinguish five or three prime utr
        for metagene in all_metagene.values():
            metagene.five_or_three_utr()
    # first last exon
    for metagene in all_metagene.values():
        metagene.first_last_exon()
    print('Done Building Metagene')
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

def pos_relative_entropy(eCLIP_prob):
    '''
    return relative entropy throughout whole transcripts
    '''
    rel_ent = {}
    for feat in featnames:
        for align in ['left', 'right']:
            pos = np.mean(np.array([eCLIP_prob['positive', feat,align, r] for r in ['rep1', 'rep2']]), axis = 0)
            neg = np.mean(np.array([eCLIP_prob['negative', feat,align, r] for r in ['rep1', 'rep2']]), axis = 0)
            rel_ent[feat, align] = density_array_entropy(pos, neg)
    return rel_ent
def remove_background(eclip, chrom, start, stop, strand, method = 'subtract'):
    '''
    remove background by comparing IP to SMInput
    return dictionary keys = rep1, rep2; values = np.array
    '''
    # how many reps are there
    rep_keys = eclip.rep_keys
    result_dict = {}

    for rep in rep_keys:
    
        ### get control density
        if 'ctrl' in eclip.read_densities.keys():
            # only 1 ctrl for both reps
            
            density_ctrl = np.nan_to_num(eclip.read_densities['ctrl'].values(chrom, start, stop, strand),0)
        else:
            # choose corresponding ctrl
            no_rep = rep.replace('rep', '')
            density_ctrl = np.nan_to_num(eclip.read_densities['ctrl'+no_rep].values(chrom, start, stop, strand),0)
        
        ### get IP density
        density_IP = np.nan_to_num(eclip.read_densities[rep].values(chrom, start, stop, strand),0)

        ### revert sign if on the negative strand
        if strand == '-':
            density_ctrl = -density_ctrl
            density_IP = -density_IP
        
        ### remove background
        if method == 'subtract':
            result = np.array(density_IP) - np.array(density_ctrl)
            # no negative value
            result[result < 0] = 0
        
        result_dict[rep] = result
    
    
    
    return result_dict
