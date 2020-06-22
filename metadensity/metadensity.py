import pyBigWig
import pandas as pd
import sys
sys.path.append('/home/hsher/rbp-maps/maps/')

from density.ReadDensity import ReadDensity
import os
from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np
basedir = '/home/hsher/seqdata/eclip_raw/'
from scipy.stats import entropy


transcript = BedTool('/home/hsher/projects/gencode_transcript.gff3')
gencode_feature = BedTool('/home/hsher/projects/gencode_combine_sorted.gff3')
featnames = ['five_utr', 'first_exon', 'exon', 'intron', 'last_exon','three_utr']

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
        
class eCLIP:
    def __init__(self, density_rep1 = None, density_rep2 = None, density_ctrl = None):
        self.rep1 = density_rep1 # Read Density Objects
        self.rep2 = density_rep2
        self.ctrl = density_ctrl

    def from_pd_series(self, series):
        ''' from the above encode data dataframe'''
        rep1, rep2, ctrl = make_density(series)
        self.__init__(rep1, rep2, ctrl)
        self.uID = series['uID'].values[0]
        self.name = series['RBP'].values[0]
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
        self.idr_transcript = genome_coord.intersect(self.idr, s=True)
    def find_negative_example(self, genome_coord = transcript):
        ''' find negative regions as no IDR neither individual peaks'''
        self.no_peak = genome_coord.intersect(self.idr, s=True, v=True).intersect(self.peak1, s = True, v =True).intersect(self.peak2, s = True, v =True)
     
    def build_metagene(self, sample_no = 200):
        ''' Use function `Build_many_metagenes()` to get regions in UTR, intron, exon for each transcripts in self.idr_transcript and self.no_peak; store metagene in self.idr_metagene/self.neg_metagene'''
        self.idr_metagene = Build_many_metagene(self.idr_transcript, sample_no = sample_no)
        self.neg_metagene = Build_many_metagene(self.no_peak, sample_no = sample_no)
    def get_metadensity(self):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object'''
        _ = [m.metadensity(self) for m in self.idr_metagene.values()]
        _ = [m.metadensity(self) for m in self.neg_metagene.values()]
    def get_quantile_metadensity(self, q = 40):
        ''' making features into quantile after averaging'''
        _ = [m.quantile_metadensity(q = q) for m in self.idr_metagene.values()]
        #_ = [m.quantile_metadensity() for m in self.neg_metagene.values()] #### negative example metadensity don't work due to not enough unique value
    def get_feature_density_array(self, feature, target_len, align, pad_value = np.nan, example = 'positive', use_quantile = False):
        ''' align features by padding np.nan to the 5' or 3' end '''
        d1 = []
        d2 = []
        if example == 'positive':
            metagenes = self.idr_metagene.values()
        else:
            metagenes = self.neg_metagene.values()
        for m in metagenes:
            
            
            
            if use_quantile:
                # see if multi-feature average is called
                if type(m.qdensities[self.uID]['rep1'][feature]) == tuple:
                    if align == 'left':
                        d1.append(trim_or_pad(m.qdensities[self.uID]['rep1'][feature][0], target_len, align, pad_value))
                        d2.append(trim_or_pad(m.qdensities[self.uID]['rep2'][feature][0], target_len, align, pad_value))
                    elif align == 'right':
                        d1.append(trim_or_pad(m.qdensities[self.uID]['rep1'][feature][1], target_len, align, pad_value))
                        d2.append(trim_or_pad(m.qdensities[self.uID]['rep2'][feature][1], target_len, align, pad_value)) 
                else: # only 1 entry, no such align problem!
                    d1.append(trim_or_pad(m.qdensities[self.uID]['rep1'][feature], target_len, align, pad_value))
                    d2.append(trim_or_pad(m.qdensities[self.uID]['rep2'][feature], target_len, align, pad_value))
            else:
                if type(m.densities[self.uID]['rep1'][feature]) == tuple:
                    if align == 'left':
                        d1.append(trim_or_pad(m.densities[self.uID]['rep1'][feature][0], target_len, align, pad_value))
                        d2.append(trim_or_pad(m.densities[self.uID]['rep2'][feature][0], target_len, align, pad_value))
                    elif align == 'right':
                        d1.append(trim_or_pad(m.densities[self.uID]['rep1'][feature][1], target_len, align, pad_value))
                        d2.append(trim_or_pad(m.densities[self.uID]['rep2'][feature][1], target_len, align, pad_value)) 
                else:
                    d1.append(trim_or_pad(m.densities[self.uID]['rep1'][feature], target_len, align, pad_value))
                    d2.append(trim_or_pad(m.densities[self.uID]['rep2'][feature], target_len, align, pad_value))
        
        # store differntly
        if use_quantile:
            self.qdensity_array[example,feature,align,'rep1']= np.stack(d1)
            self.qdensity_array[example,feature,align,'rep2']= np.stack(d2)
        else:
            self.density_array[example,feature,align,'rep1']= np.stack(d1)
            self.density_array[example,feature,align,'rep2']= np.stack(d2)
                
        
    def get_density_array(self, five_utr_len=100, three_utr_len=150, intron_len = 1500, exon_len = 150, use_quantile = False):
        ''' extract metadensity from each metagene, zero pad or trim to get np array '''
        if use_quantile:
            examples = ['positive'] # quantile can only have positive arrays
        else:
            examples = ['positive', 'negative']
        
        self.density_array = {}
        self.qdensity_array = {}
        for feature, l in zip(featnames, [five_utr_len, exon_len,exon_len, intron_len, exon_len,three_utr_len]):
            for example in examples:
                for align in ['left', 'right']:
                    #print(feature, example, align)
                    self.get_feature_density_array(feature, l, align, example = example, use_quantile = use_quantile)
    
    
    def concat_density_array(self, example = 'positive', rep = 'rep1', quantile = False):
        ''' return concatenated density array by sequence of featnames '''
        if quantile:
            concat = np.concatenate([self.qdensity_array[example, feat, align, rep] 
                        for feat in featnames for align in ['left', 'right']],axis = 1)
        else:
            
            concat = np.concatenate([self.density_array[example, feat, align, rep] 
                        for feat in featnames for align in ['left', 'right']],axis = 1)
        return concat
    def scale_density_array(self, method = 'norm_by_sum', quantile = False):
        ''' scale density array by only the features considered '''
        self.scaled_density_array = {}
        if quantile:
            examples = ['positive']
            denarray = self.qdensity_array
        else:
            examples = ['positive', 'negative']
            denarray = self.density_array
            
        for example in examples:
        
            rep1_concat = self.concat_density_array(rep = 'rep1', example = example, quantile = quantile)
            rep2_concat = self.concat_density_array(rep = 'rep2', example = example, quantile = quantile)
        
            if method == 'norm_by_sum':
                scale1 = np.nansum(rep1_concat, axis = 1)[:, None]
                scale2 = np.nansum(rep2_concat, axis = 1)[:, None]
            
            if method == 'max_scalar':
                scale1 = np.nanmax(rep1_concat, axis = 1)[:, None]
                scale2 = np.nanmax(rep2_concat, axis = 1)[:, None]
            for rep, scale in zip(['rep1', 'rep2'], [scale1, scale2]):
            
                for align in ['left', 'right']:
                    for feature in featnames:
                        self.scaled_density_array[example,feature,align,rep] = denarray[example,feature,align,rep]/scale
            
           
        
    def RBP_centric_approach(self, series, sample_no = 200):
        ''' create eCLIP object, get peaks, get examples, metagene and metadensity with one function. return eCLIP object'''
        self.from_pd_series(series)
        print('adding peaks')
        self.add_peaks()
        print('finding negative/positive examples')
        self.find_idr_transcript()
        self.find_negative_example()
        print('Building metagene and metadensity')
        self.build_metagene(sample_no = sample_no)
        self.get_metadensity()
        
        

########################## Metagene class #####################################

class Metagene:
    def __init__(self, esnt, chro, start, end, strand):
        self.ensembl_id = esnt
        self.chrom = chro
        self.start = start
        self.stop = end
        self.strand = strand
        self.five_utr = set()
        self.three_utr = set()
        self.intron = set()
        self.exon = set()
        self.densities = {} # key: eCLIP uID, values: metadensity
        self.qdensities  = {}
    def multi_feature_avg(self, feature, eCLIP, align = 'left', max_len = None, quantile = False, bins1 = None, bins2 = None):
        '''
        average and align (zero padding) multiple intron/exon
        return nanmean for both replicates
        '''
        den1 = []
        den2 = []
        
        # feature length
        if max_len == None:
            max_len = max([f[1]-f[0] for f in feature])
        
        for f in feature:
            minus1, minus2 = subtract_input(eCLIP, self.chrom, f[0], f[1], self.strand)
            
            if quantile:
                minus1 = np.digitize(minus1, bins1)
                minus2 = np.digitize(minus2, bins2)
                
            
        
            den1.append(trim_or_pad(minus1, max_len, pad_value = np.nan, align = align))
            den2.append(trim_or_pad(minus2, max_len, pad_value = np.nan, align = align))
        return np.nanmean(np.stack(den1), axis = 0), np.nanmean(np.stack(den2), axis = 0)
        
    def first_last_exon(self):
        ''' find first, last exon '''
        # seperate first and last exon
        min_start = min([e[0] for e in list(self.exon)])
        max_start = max([e[1] for e in list(self.exon)])
        
        self.first_exon = set([e for e in list(self.exon) if e[0] == min_start])
        self.last_exon = set([e for e in list(self.exon) if e[1] == max_start])
        
        self.exon = self.exon - self.first_exon - self.last_exon
        
    def concat_density(self, uID, rep):
        ''' since some feature would have left/right, we need a special function to return concat density '''
        left_densities = [feat[0] if type(feat) == tuple else feat for feat in self.densities[uID][rep].values()]
        return np.concatenate(left_densities)
        
        
    def metadensity(self, eCLIP):
        ''' get metadensity from eCLIP object (containing density) 
        store. self.densities[eCLIP.uid] as dictionary {'rep1': exon: np.array}'''
                    
        # call first and last exon
        self.first_last_exon()       
        
        # initialize
        self.densities[eCLIP.uID] = {}
        self.densities[eCLIP.uID]['rep1'] = {}
        self.densities[eCLIP.uID]['rep2'] = {}
                
        # get average density per feature
        for feature, fname in zip([self.five_utr, self.first_exon, self.exon, self.intron, self.last_exon, self.three_utr], featnames):
            if len(feature) == 0: # no such feature
                minus1 = np.empty(1) # np.empty does not always yield nan
                minus2 = np.empty(1)
                minus1[:] = np.nan
                minus2[:] = np.nan
                
            elif len(feature) == 1:
                feature = list(feature)[0]
                minus1, minus2 = subtract_input(eCLIP, self.chrom, feature[0], feature[1], self.strand)
                minus1 = np.array(minus1)
                minus2 = np.array(minus2)
            else:
                minus1_l, minus2_l = self.multi_feature_avg(feature, eCLIP, align = 'left') ######## automatically align to left, causing low value in distal intron?
                minus1_r, minus2_r = self.multi_feature_avg(feature, eCLIP, align = 'right')
                minus1 = (minus1_l, minus1_r)
                minus2 = (minus2_l, minus2_r)
            
            self.densities[eCLIP.uID]['rep1'][fname] = minus1
            self.densities[eCLIP.uID]['rep2'][fname] = minus2
        
        # normalization
        denom1 = np.nansum(self.concat_density(eCLIP.uID, 'rep1'))
        denom2 = np.nansum(self.concat_density(eCLIP.uID, 'rep2'))
        #print(denom1, denom2) # debug
        
        if denom1 != np.nan and denom1 != 0:
            for fname in featnames:
                feat = self.densities[eCLIP.uID]['rep1'][fname]
                if type(feat) == tuple:
                    self.densities[eCLIP.uID]['rep1'][fname] = (feat[0]/denom1, feat[1]/denom1)
                else:
                    self.densities[eCLIP.uID]['rep1'][fname] = feat/denom1
        
        if denom2 != np.nan and denom2 != 0:
            for fname in featnames:
                feat = self.densities[eCLIP.uID]['rep2'][fname]
                if type(feat) == tuple:
                    self.densities[eCLIP.uID]['rep2'][fname] = (feat[0]/denom2, feat[1]/denom2)
                else:
                    self.densities[eCLIP.uID]['rep2'][fname] = feat/denom2
        
    
        
        
    def quantile_metadensity(self, q = 40):
        
        ''' quantile metadensity to equal sized bins at transcript level'''
        
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
        ''' quantilize metadensity uning raw density (not pooled, not averaged)'''
        uID = eCLIP.uID
        
        # initialize empty dictionary
        self.qdensities[uID] = {}
        feature_den1 = []
        feature_den2 = []
        
        # get bins
        transcript1, transcript2 = subtract_input(eCLIP, self.chrom, self.start, self.stop, self.strand)
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
                minus1, minus2 = subtract_input(eCLIP, self.chrom, feature[0], feature[1], self.strand)
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
    ''' Create Metagene object for regions for key transcript '''
    if len(key_transcript)< sample_no:
        n = len(key_transcript)
    else: 
        n = sample_no
    
    # put into transcript
    all_metagene = {}
    for i in range(n):
        enst = key_transcript[i].fields[-1].split(';')[3].replace('transcript_id=', '')
        metagene = Metagene(enst, key_transcript[i].chrom, key_transcript[i].start, key_transcript[i].stop, key_transcript[i].strand)
        all_metagene[enst] = metagene
    
    # write where is intron, exon...
    for i in gencode_feature.filter(lambda x: x.fields[-1].split(';')[3].replace('transcript_id=', '') in all_metagene.keys()):
        enst = i.fields[-1].split(';')[3].replace('transcript_id=', '')
        feature_type = i.fields[2]
        if feature_type == 'transcript':
            all_metagene[enst].intron.update([(i.start, i.stop)])
        if feature_type == 'exon':
            all_metagene[enst].exon.update([(i.start, i.stop)])
        if feature_type == 'five_prime_UTR':
            all_metagene[enst].five_utr.update([(i.start, i.stop)])
        if feature_type == 'three_prime_UTR':
            all_metagene[enst].three_utr.update([(i.start, i.stop)])
    return all_metagene

########################################## calculate meta-density ##########################################
def subtract_input(eclip, chrom, start, stop, strand):
    '''
    return IP density minus input density for 2 biological replicates
    '''
    density1 = np.nan_to_num(eclip.rep1.values(chrom, start, stop, strand),0)
    density2 = np.nan_to_num(eclip.rep2.values(chrom, start, stop, strand),0)
    
    
    
    # get input
    density_ctrl = np.nan_to_num(eclip.ctrl.values(chrom, start, stop, strand),0)
    
    if strand == '-':
        density1 = -density1
        density2 = -density2
        density_ctrl = -density_ctrl
    
    # transcript level normalization of density (if nansum = 0, that means there is no signal)(NOT CORRECT WHEN THIS IS CALLED ON A FEATURE LEVEL)
    #if np.nansum(np.array(density1)) != 0:
    #    norm_den1 = np.array(density1)/np.nansum(np.array(density1))
    #else:
    #    norm_den1 = np.array(density1)
    #if np.nansum(np.array(density2)) != 0:
    #    norm_den2 = np.array(density2)/np.nansum(np.array(density2))
    #else:
    #    norm_den2 = np.array(density2)
    #if np.nansum(np.array(density_ctrl)) != 0:
    #    norm_ctrl = np.array(density_ctrl)/np.nansum(np.array(density_ctrl))
    #else:
    #    norm_ctrl = np.array(density_ctrl)
    
    # substract after normalization
    #minus1 = norm_den1 - norm_ctrl
    #minus2 = norm_den2 - norm_ctrl
    
    minus1 = np.array(density1) - np.array(density_ctrl)
    minus2 = np.array(density2) - np.array(density_ctrl)
    
    # no negative value 
    minus1[minus1 < 0] = 0
    minus2[minus2 < 0] = 0
    
    
    return minus1, minus2

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
