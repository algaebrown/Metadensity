
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
import deepdish as dd


from metadensity.readdensity import ReadDensity
from metadensity.truncation import read_start_sites
# from . import settings
from metadensity.config import settings
transcript = BedTool(settings.transcript_fname) #TODO
gencode_feature = BedTool(settings.gencode_feature_fname)
print('Using: ', settings.transcript_fname)


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
    def __init__(self, name='RBP', uID=None, single_end=True, rep_keys = ['rep1', 'rep2'], read_densities={}, peaks = {}, idr = None, read2 = True):
        # experimental variables
        self.name = name
        if uID:
            self.uID = uID
        else:
            self.uID = name
        self.rep_keys = rep_keys
        self.single_end = single_end
        self.read2 = read2

        # data
        self.read_densities = read_densities # bigwig, bam
        self.peaks = peaks # bed 
        self.idr = idr # bed
    
    @classmethod
    def from_series(cls, series, single_end = True, read2 = True):
        """[build eCLIP object from series containing absolute path for .bam, .bw, .beds,
        ]

        Args:
            series ([pd.Series]): [RBP for name; uid for uID.
        bam_0, bam_1, bam_control (or bam_control_0, bam_control_1);
        plus_0, plus_1, plus_control for pos.bw; minus_0, minus_1, minus_control for neg.bw;
        bed_0, bed_1 for individual peaks, idr for IDR peaks.
        _0, _1 is for replicate 1,2; you can extend to the # of reps you have.
        for control, if you have more than 1 SMInput, use bam_control_0, bam_control_1...etc.
        If only 1 SMInput, simply omit _0. use bam_control, plue_control...etc]
            single_end (bool, optional): [if your CLIP is done in single end set to True else Metatruncate will have problem]. Defaults to True.
            read2 (bool, optional): [Use read 2 only if sequencing is paired-end.]. Defaults to True.

        Returns:
            [eCLIP]: [a object that holds eCLIP data]
        """
        ''' :
        
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
        if 'idr' in series.index:
            idr = BedTool(series['idr'])
        else:
            idr = None

        return cls(name=series['RBP'], uID=series['uid'], single_end=single_end, rep_keys = rep_keys, read_densities=read_densities, peaks = peaks, idr = idr, read2 = read2)

        

    
    def find_idr_transcript(self, genome_coord = transcript):
        
        """[find transcript containing at least 1 IDR peak]
        Kwargs:
            genome_coord: (BedTool) transcript coordintae. default are canonical transcripts.
        """
        if self.idr:
            self.idr_transcript = genome_coord.intersect(self.idr, s=True, u = True).saveas()
        else:
            self.idr_transcript = genome_coord.intersect(self.peaks[self.rep_keys[0]], s=True, u = True).saveas()
            print('Warning: No IDR file input, falling back to using peaks from rep {}'.format(self.rep_keys[0]))

    ##################### FOR ENCODE DATA USE ONLY ###################################
    def enough_transcripts(self, coord = transcript, thres = 200, key = 'rep1', attr_to_return = 'ID'):
        """[Filter for transcript id that has enough reads in IP]

        Args:
            coord ([BedTool], optional): [The intervals to count number of reads. For example, if set to transcript, 
            will return transcript id with at least "thres" reads. Has to be in gff3. Check if interval.attrs['ID'] exist. ]. Defaults to transcript.
            thres (int, optional): [the number of reads]. Defaults to 200.
            key (str, optional): [library to look at. Set to any keys in self.read_densities]. Defaults to rep1.

        Returns:
            [list]: [list of interval ids]
        """
        enough_read = set()
        for t in coord:
            if len(list(self.read_densities[key].bam.fetch(t.chrom, t.start,t.end))) > thres:
                enough_read.add(t.attrs[attr_to_return])
        return enough_read
    def enough_coverage_transcript(self, coord = transcript, thres = 0.001, key = 'rep1', attr_to_return = 'ID'):
        """[Filter for transcript id that has enough coverage in IP library]

        Args:
            coord ([BedTool], optional): [The intervals to count number of reads. For example, if set to transcript, 
            will return transcript id with at least "thres" reads. Has to be in gff3. Check if interval.attrs['ID'] exist. ]. Defaults to transcript.
            thres (int, optional): [the number of reads]. Defaults to 200.
            key (str, optional): [library to look at. Set to any keys in self.read_densities]. Defaults to rep1.

        Returns:
            [list]: [list of interval ids]
        """
        enough_read = set()
        for t in coord:
            coverage=len(list(self.read_densities[key].bam.fetch(t.chrom, t.start,t.end)))/(t.end-t.start)
            if coverage > thres:
                enough_read.add(t.attrs[attr_to_return])
        return enough_read
    def calculate_coverage(self, coord = transcript):
        """[Calculate coverage for all bams for intervals in coord]
        Args:
            coord ([BedTool], optional): [The intervals to count number of reads. For example, if set to transcript, 
        Returns:
            coveragedf
        """
        bam_fnames = []
        libs = list(self.read_densities.keys())
        for key in libs:
            bam_fnames.append(self.read_densities[key].bam.filename.decode('ascii')) #byte from pysam
            
        read_count = coord.multi_bam_coverage(bams=bam_fnames, s = True).to_dataframe()
        columns = list(read_count.columns)

        
        columns[-len(libs):]= libs
        read_count.columns = columns

        read_count['length'] = read_count['end']-read_count['start']

        read_count[[f'{l}_rpk' for l in libs]]=read_count[libs].div(read_count['length'], axis = 0)*1000
        # calculate read per bp
        return read_count
### Metadensity

class Meta:
    ''' superclass of to inherit from '''
    def __init__(self, eCLIP, name, sample_no = 200, metagenes = None, transcripts = None, transcript_ids = None, deep_dish_path = None):
        
        self.eCLIP = eCLIP
        self.name = name

        if deep_dish_path:
            # load from precomputed data!
            self.density_array = dd.io.load(deep_dish_path)
            self.featnames = list(set([k[0] for k in self.density_array.keys()]))
        

        else:
            # compupte from scratch
            if metagenes:
                self.metagene = metagenes # precomupted 
            elif transcript_ids:
                self.build_metagene(sample_no = sample_no, tids = transcript_ids)
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
            
            self.density_array = {}

    
    
    def build_metagene(self, sample_no = 200, tids = None):
        ''' Use function `Build_many_metagenes()` to get regions in UTR, intron, exon for each transcripts in self.idr_transcript and self.no_peak; store metagene in self.idr_metagene/self.neg_metagene'''
        if not tids:
            tids = [s.attrs['transcript_id'] for s in self.transcript]
        self.metagene = Build_many_metagene(tids, sample_no = sample_no)
    def get_feature_density_array(self, feature, target_len, align, pad_value = np.nan, use_truncation = False):
        """[align features by padding np.nan to the 5' or 3' end]

        Args:
            feature ([str]): [type of feature to fetch. keys of metagene.features, can be 'exon', 'first_exon', 'intron'...etc. check metagene.features]
            target_len ([int]): [The windowing length of that feature.]
            align ([str]): [whether to align to 5 prime or 3 prime, 'left'(5 prime) or 'right'(3 prime)]
            pad_value ([float], optional): [What value to pad if the feature is shorter than designated target_len]. Defaults to np.nan.
            use_truncation (bool, optional): [Whether to use truncation or not. If use truncation, fetch values from metagene.truncations. else fetch from metagene.densities]. Defaults to False.
        """
        
        metagenes = self.metagene.values()
        for rep in self.eCLIP.rep_keys:
            rep_densities = []
            for m in metagenes:          
                        
                
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
            self.density_array[feature,align,rep]= np.stack(rep_densities)
    
    
    def get_density_array(self, use_truncation = False, full_CDS = True):
        """[Extract density for all features in all metagene, wrap around self.get_feature_density_array]

        Args:
            use_truncation (bool, optional): [Whether to use truncation or not. If use truncation, fetch values from metagene.truncations. else fetch from metagene.densities]. Defaults to False.
            full_CDS (bool, optional): [Whether to concat all CDS and make a feature called 'full_CDS']. Defaults to True.
        """
        if full_CDS:
            # build CDS feature
            try:
                _ = [metag.build_full_cds(self.eCLIP, truncate = use_truncation) for metag in self.metagene.values()]
            except Exception as e:
                print(e)

        for feature in self.featnames:
            for align in ['left', 'right']:
                flen = settings.feat_len[feature]
                self.get_feature_density_array(feature, flen , align, use_truncation = use_truncation)
    def concat_density_array(self, rep = 'rep1', features = ['exon', 'intron']):
        """[concatenate density array along specified features. For examples, exon is (n * 200), intron is (n * 500), then the concat array will be (n * 700)]

        Args:
            rep (str, optional): [which replicate to concat]. Defaults to 'rep1'.
            features (list, optional): [what features to concat]. Defaults to 'self.featnames'.
            
        Returns:
            [np.array]: [concatenated array]
        """
        ''' return concatenated density array by sequence of featnames '''
        concat = np.concatenate([self.density_array[feat, align, rep] 
                    for feat in features for align in ['left', 'right']]
                    ,axis = 1)
        return concat
    def scale_density_array(self, method = 'norm_by_sum'):
        """[scale density array by only the features considered]

        Args:
            method (str, optional): ['norm_by_sum', 'max_scalar']. Defaults to 'norm_by_sum'.
        
        Returns:
            [dict]: [dictionary, keys = [feat, align, rep], values = np.array]
        """
        
        scaled_density_array = {}
        for rep in self.eCLIP.rep_keys:
            rep_concat = self.concat_density_array(rep = rep)
            if method == 'norm_by_sum':
                scale = np.nansum(rep_concat, axis = 1)[:, None]
            if method == 'max_scalar':
                scale = np.nanmax(rep_concat, axis = 1)[:, None]
            for align in ['left', 'right']:
                for feature in self.featnames:
                    scaled_density_array[feature,align,rep] = self.density_array[feature,align,rep]/scale
        return scaled_density_array
    def apply(self, func, axis = 0):
        """[apply function to all features in self.density_array ]

        Args:
            func ([fuction]): [examples likenp.nanmean, np.median, np.std]
            axis (int, optional): [along with axis do the function applies to]. Defaults to 0.

        Returns:
            [dict]: [dictionary, keys = [feat, align, rep], values = np.array]
        """
        
        result = {}
        for key in self.density_array.keys():
            result[key] = func(self.density_array[key], axis = axis)
        return result
    
    def convolve(self, filter, axis = 0):
        """[wrapper around np.convolve, to run in all features,  use to smooth out signals ]

        Args:
            filter ([np.array]): [v in np.convolve]
            axis (int, optional): [which dimension to do]. Defaults to 1.

        Returns:
            [dict]: [dictionary, keys = [feat, align, rep], values = np.array]
        """
        ''' wrapper around np.convolve,  use to smooth out signals '''
        result = {}
        for key in self.density_array.keys():
            result[key] = np.apply_along_axis(lambda m: np.convolve(m, filter, mode='full'), axis=axis, arr = self.density_array[key])
        return result
    def save_deepdish(self, outdir):
        """[save self.density_array to deepdish]

        Args:
            outdir ([str]): [path to save]
        """
        dd.io.save(outdir, self.density_array)



class Metatruncate(Meta):
    """[Metagene for 5' read start]

    Args:
        Meta ([cls]): [inherit from superclass Meta]
    """
    def __init__(self, eCLIP, name, sample_no = 200, background_method = 'subtract', normalize = True, metagenes = None, transcripts = None, transcript_ids = None, deep_dish_path = None):
        """[Create a Metatruncate class (using 5 prime read truncation)]

        Args:
            eCLIP ([eCLIP]): [eCLIP object containing the .bam .bw]
            name ([str]): [The name of this object. Will be used for plotting]
            sample_no (int, optional): [Number of transcripts to calculate. If is shorter than the transcripts specifies, then it will use the first n]. Defaults to 200.
            background_method (str, optional): ['subtract'; 'subtract normal'; 'relative information']. Defaults to 'subtract'.
            normalize (bool, optional): [True to use % of edits in position; False to use raw background subtract signals]. Defaults to True.
            metagenes ([dictionary of metagenes], optional): [Use pre-built metagene. Will override sample_no, transcripts and transcript_ids]. Defaults to None.
            transcripts ([BedTool], optional): [Use BedTool containing transcripts]. Defaults to None.
            transcript_ids ([list], optional): [list of transcript or gene ids. the ids needs to be in the pre-built coodinate]. Defaults to None.
            deep_dish_path ([str], optional): [path to precomputed density array deepdish.h5. If specified, directly load from path and override all options]. Defaults to None.
        """
        # generate metagene coords
        super().__init__(eCLIP, name, sample_no, metagenes = metagenes, transcripts = transcripts, transcript_ids = transcript_ids,
        deep_dish_path = deep_dish_path)

        # automatically run
        if not deep_dish_path:
            self.get_truncation(background_method = background_method, normalize = normalize)
        self.background_method = background_method
        self.normalize = normalize
    
    
    def get_truncation(self, background_method = 'subtract', normalize = True):
        """[calculate metadensity for each metagene.
        store density in each metagene object]

        Args:
            background_method (str, optional): ['subtract'; 'subtract normal'; 'relative information']. Defaults to 'subtract'.
            normalize (bool, optional): [ True to use % of edits in position; False to use raw background subtract signals]. Defaults to True.
        """
        _ = [m.get_average_feature(self.eCLIP, background_method = background_method, normalize = normalize, truncate = True) for m in self.metagene.values()]
    
    
    
class Metadensity(Meta):
    """[Metagene using read coverage (RBP protected fragments)]

    Args:
        Meta ([cls]): [superclass]
    """
    def __init__(self, eCLIP, name, sample_no = 200, background_method = 'subtract', normalize = True, metagenes = None, transcripts = None, transcript_ids = None, deep_dish_path = None):
        """[Create a Metadensity class (using read coverage)]

        Args:
            eCLIP ([eCLIP]): [eCLIP object containing the .bam .bw]
            name ([str]): [The name of this object. Will be used for plotting]
            sample_no (int, optional): [Number of transcripts to calculate. If is shorter than the transcripts specifies, then it will use the first n]. Defaults to 200.
            background_method (str, optional): ['subtract'; 'subtract normal'; 'relative information']. Defaults to 'subtract'.
            normalize (bool, optional): [True to use % of edits in position; False to use raw background subtract signals]. Defaults to True.
            metagenes ([dictionary of metagenes], optional): [Use pre-built metagene. Will override sample_no, transcripts and transcript_ids]. Defaults to None.
            transcripts ([BedTool], optional): [Use BedTool containing transcripts]. Defaults to None.
            transcript_ids ([list], optional): [list of transcript or gene ids. the ids needs to be in the pre-built coodinate]. Defaults to None.
            deep_dish_path ([str], optional): [path to precomputed density array deepdish.h5. If specified, directly load from path and override all options]. Defaults to None.
        """
        # generate metagene coords
        super().__init__(eCLIP, name, sample_no, metagenes = metagenes, transcripts = transcripts, transcript_ids = transcript_ids
        ,deep_dish_path = deep_dish_path)

        # automatically run
        if not deep_dish_path:
            self.get_metadensity(background_method = background_method, normalize = normalize)
        self.background_method = background_method
        self.normalize = normalize
        
        
    def get_metadensity(self, background_method = 'subtract', normalize = True):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object'''
        _ = [m.get_average_feature(self.eCLIP, background_method = background_method, normalize = normalize, truncate = False) for m in self.metagene.values()]
    
    

    

        
        

########################## Metagene class #####################################
class Window:
    def __init__(self, start, end):
        self.start = start
        self.end = end

        self.start_name = None # set for plotting
        self.end_name = None
class Metagene:
    """[A class to hold a gene and its feature]
    """
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        """[Create a Metagene object]

        Args:
            esnt ([string]): [ID of the gene. store in Metagene.id]
            chro ([string]): [chromosome. store in Metagene.chrom]
            start ([int]): [start of the gene/transcript. self.start]
            end ([int]): [end of the gene/transcript. self.stop]
            strand ([string]): [has to be '+' or '-']
            transcript_type ([string]): [type of transcript]
            features ([list], optional): [description]. Defaults to None.
        Attributes:
            self.featnames ([list]): [for things to appear in density array, the name of the feature needs to be this list.]
            self.coverage ([dict]): [raw coverage values from bigwig]
            self.sites ([dict]): [raw truncation values from bam]
            self.value([dict]): [background removed, normalized value]
            self.densities ([dict])= [windowing for each feature from self.value derived from coverage]
            self.truncations ([dict])= [windowing for each feature from self.value derived from truncation]
        """
        self.id = esnt # can be other depending on stuffs
        self.type = transcript_type
        self.chrom = chro
        self.start = start
        self.stop = end
        self.strand = strand

        if not features:
            self.features = defaultdict(set)
        else:
            self.features = defaultdict(set, features)
        
        ## These needs to be downloaded
        self.featnames = ['first_exon', 'exon', 'intron', 'branchpoint', 'branchpoint_pred','last_exon', 'polyAsite', 'polyAsignal']
        
        

        # TODO record raw value before pooling, save I/O time
        self.coverage = {} # key: eCLIP uID, values: coverage RAW from 5' end to 3' end
        self.sites = {} # key: eCLIP uID, values: RAW truncation count
        
        # normalized of entire transcript
        self.value = {}

        # pooled metagene density from such
        self.densities = {} # key: eCLIP uID, values: metadensity
        self.truncations = {}
    
    @classmethod
    def from_dict(cls, coord_dict):
        """[class method, generate metagene object from precomputed coordinates]

        Args:
            coord_dict ([type]): [description]

        Returns:
            [type]: [description]
        """
        return cls(coord_dict['id'],
        coord_dict['chrom'], coord_dict['start'], coord_dict['end'], coord_dict['strand'],
         coord_dict['type'], features = coord_dict['features'])
    ################################## about annotation ###########################################
    
    
    def order_multi_feature(self, feature = 'exon'):
        """[order multifeature from five prime to three prime end]

        Args:
            feature (str, optional): [type of feature to order. 'exon' or 'cds']. Defaults to 'exon'.

        Returns:
            [list]: [a list of features, 5 prime to 3 prime]
        """
        
        feat_list = list(self.features[feature])
        if self.strand == '-':
            # the largest 5 prime end is the first
            
            sorted_middle_exon = sorted(feat_list, key = lambda interval: interval[1], reverse = True)
        else:
            sorted_middle_exon = sorted(feat_list, key = lambda interval: interval[0], reverse = False)
        return list(self.features['first_'+feature])+ sorted_middle_exon + list(self.features['last_'+feature])
    
    def concat_multi_feature(self, value, feature = 'exon'):
        """[concatentae multiple feature in order, for example, multiple exon/CDS, from 5' to 3']

        Args:
            value ([np.array]): [the original value to slice from]
            feature (str, optional): [type of feature, 'exon' or 'cds']. Defaults to 'exon'.

        Returns:
            [np.array]: [concatenate density of multiple feature]
        """
        feat_in_order = self.order_multi_feature(feature = feature)
        concat_val = []
        for feat in feat_in_order:
            relative = self.to_relative_axis(feat)
            if self.strand == '-':
                concat_val.append(value[relative[0]: relative[1]][::-1]) # TODO check if we need to +1
            else:
                concat_val.append(value[relative[0]: relative[1]]) # don't need to +1
        return np.concatenate(concat_val)
    
    def concat_density(self, uID, rep, truncation = False):
        ''' since some feature would have left/right, we need a special function to return concat density '''
        if truncation:
            left_densities = [feat[0] if type(feat) == tuple else feat for feat in self.truncations[uID][rep].values()]
        else:
            
            left_densities = [feat[0] if type(feat) == tuple else feat for feat in self.densities[uID][rep].values()]
        return np.concatenate(left_densities)

    def create_feature(self, interval, feature_name, length=50):
        """[create new feature in metagene]

        Args:
            interval ([int or tuple]): [if int, create point feature at position x; if tuple of (int, int), create interval.Absolute value in genome coordinate]
            feature_name ([str]): [name of the feature]
        Example:
            branchpoint for this gene is 500. The gene is from (100, 300).
            create_feature(interval=500, feature_name = 'branchpoint')

            a snoRNA is nested in this gene from (130, 150).
            call create_feature(interval(130,150), feature_name = 'snoRNA')
        """
        try:
            self.features[feature_name].add(
                (
                int(interval[0]), 
                int(interval[1])
                )
                )
        except:
            # CREATING point feature
            self.features[feature_name].add(
                int(interval)
                )
        self.featnames.append(feature_name)
        settings.feat_len[feature_name] = length
    def create_downstream_feature(self, query_interval, feature_type = 'intron', feature_name = 'terminal_exon'):
        """[create closest, 3' feature of the query interval
        for example, intron list: (100,200), (400,500), (700,900), query exon (500,700), return (700,900) if '+' else (400,500)]

        Args:
            query_interval ([tuple]): [(int,int), the interval if interest]
            feature_type (str, optional): ['exon','intron']. Defaults to 'intron'.
            feature_name (str, optional): [name of the query feature]. Defaults to 'terminal_exon'.

        Returns:
            [None]: [nothing]
        """

        orders = self.order_multi_feature(feature_type) # order five prime to three prime
        for o in orders:
            if self.strand == '+':
                if o[0]>=query_interval[1]: # end of the query is 5' to start of feature
                    next_feat = o
                    break
            else:
                if o[1]<= query_interval[0]:
                    next_feat = o
                    break
        try:
            self.create_feature(next_feat, f'{feature_name}_downstream_{feature_type}')
            
        except:
            self.featnames.append(f'{feature_name}_downstream_{feature_type}')
    def create_upstream_feature(self, query_interval, feature_type = 'intron', feature_name = 'terminal_exon'):
        """[create closest, 3' feature of the query interval
        for example, intron list: (100,200), (400,500), (700,900), query exon (500,700), return (700,900) if '+' else (400,500)]

        Args:
            query_interval ([type]): [description]
            feature_type (str, optional): [description]. Defaults to 'intron'.
            feature_name (str, optional): [description]. Defaults to 'terminal_exon'.

        Returns:
            [type]: [description]
        """
        

        orders = self.order_multi_feature(feature_type)[::-1] # order three prime to five prime
        for o in orders:
            if self.strand == '+':
                if o[1]<=query_interval[0]: # end of the query is 5' to start of feature
                    next_feat = o
                    break
            else:
                if o[0]>= query_interval[1]:
                    next_feat = o
                    break
        try:
            self.create_feature(next_feat, f'{feature_name}_upstream_{feature_type}')
            
        except:
            self.featnames.append(f'{feature_name}_upstream_{feature_type}')
    
    def to_relative_axis(self, interval):
        """[convert absolute genome position to relative axis, 0=five prime end of transcript/gene. Ex: self.start = 1000, self.stop = 5000, self.strand = -, interval = (1000,1005)
        return (5000-1005, 5000-1000) = (3995, 4000)
        if strand +, (1000-1000, 1005-1000)]

        Args:
            interval ([tuple]): [tuple of (int, int), absolute genome position]

        Returns:
            [tuple]: [relative position of that interval, tuple of (int, int)]
        """
        
        if self.strand == '-':
            relative_feature = (self.stop - interval[1], self.stop - interval[0])
        else:
            relative_feature = (interval[0]-self.start, interval[1] - self.start)
        return relative_feature

    ################################## about raw values ###########################################
    def truncation_count(self, eCLIP, rep, interval):
        """[fetch 5 prime end of read from the bam file. return a list of read start count for each position in feature
        ex: sites return [10,11,12,13]; feature = (10,15); return [1,1,1,1,0] count from 5' to 3']

        Args:
            eCLIP ([eCLIP]): [eCLIP object]
            rep ([str]): ['rep1']
            interval ([tuple]): [(int, int), absolute genome position of things]

        Returns:
            [np.array]: [read start for each position of designated interval, 5 prime to 3 prime]
        """
        
        if isinstance(eCLIP, STAMP):
            s = '{}\t{}\t{}\t{}\t{}\t{}'.format(self.chrom, str(interval[0]), str(interval[1]), '.', '.', self.strand)
            bed = BedTool(s, from_string=True).saveas()
            sites_bed = eCLIP.peaks[rep].intersect(bed, s = True).saveas()
            sites = [int(s.start) for s in sites_bed]
            
        else:
            
            sites = read_start_sites(eCLIP.read_densities[rep].bam, chrom=self.chrom, start=interval[0], end=interval[1], strand=self.strand, single_end = eCLIP.single_end, read2 = eCLIP.read2)
        
        site_count = Counter(sites)
        
        pos_count = []
        if self.strand == '+':
            for pos in range(interval[0], interval[1]+1):
                pos_count.append(site_count[pos])
        if self.strand == '-':
            for pos in range(interval[1], interval[0]-1, -1): # 15, 14, 13, 12, 11, 10
                pos_count.append(site_count[pos])
        return np.array(pos_count)  
    
    def coverage_value(self, eCLIP, rep, interval):
        """[return raw RPM value from bigwig files for each position in feature 
        from 5' to 3' end]

        Args:
            eCLIP ([eCLIP]): [eCLIP object]
            rep ([str]): ['rep1']
            interval ([tuple]): [(int, int), absolute genome position of interval]

        Returns:
            [np.array]: [values from bigwig, 5 prime to 3 prime]
        """
        
        rpm = np.nan_to_num(eCLIP.read_densities[rep].values(self.chrom, interval[0], interval[1], self.strand),0)
        # the Read density object takes care of strandedness
        if self.strand == '-' and np.any(rpm<0): ## default the values are negative for those from the eCLIP bigwigs
            #TODO
            rpm = -rpm
        
        return rpm
    
    def get_raw_value(self, eCLIP, truncate = False):
        """[fetch raw coverage/truncation value for whole gene 
        write in self.coverage or self.sites]

        Args:
            eCLIP ([eCLIP]): [eCLIP value]
            truncate (bool, optional): [use truncation or not. if not, fetch from bigwig]. Defaults to False.

        Returns:
            [type]: [description]
        """
        
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
        """[remove background by method
        method = 'subtract': subtract IP from Input
        method = 'subtract normal': subtract IP distribution from Input distribution
        method = 'relative information': take relative entropy
        method = 'fold change': log2(IP/IN)
        method = 'get null': get input density]

        Args:
            eCLIP ([eCLIP]): [description]
            rep ([str]): [description]
            method (str, optional): [description]. Defaults to 'subtract'.
            truncate (bool, optional): [description]. Defaults to False.

        Returns:
            [type]: [description]
        """
        
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
            
        if method == 'relative information' or method == 'fold change':
            # handle 0, add pseudocount
            
            if truncate:
                ip_pseudocount = 0.0001 #TODO this is not compatible with density
                control_pseudocount = 0.0001
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
            
            

            if method == 'relative information':
                values = np.multiply(
                    ip_norm,
                    np.log2(np.divide(ip_norm, input_norm))
                    ) # pi * log2(pi/qi)
            elif method == 'fold change':
                values = np.log2(np.divide(ip_norm, input_norm))
        if method == None:
            values = raw # do nothing with background
        
        return values   
    
    
    ################################## about annotation ###########################################
    
    

        
    
        

    def get_value(self, eCLIP, background_method = 'subtract', normalize = True, truncate = True):
        """[perform background removal, normalization and store in self.value[eCLIP.uID][rep_key]]

        Args:
            eCLIP ([type]): [description]
            background_method (str, optional): [description]. Defaults to 'subtract'.
            normalize (bool, optional): [description]. Defaults to True.
            truncate (bool, optional): [description]. Defaults to True.
        """
        
        # fetch raw sites for all reps and ctrl
        self.get_raw_value(eCLIP, truncate = truncate)

        # different raw signals are used
        if truncate:
            raw = self.sites
        else:
            raw = self.coverage

        # initialize emtpy dictionary
        self.value[eCLIP.uID] = {}

        # deal with background
        for rep in eCLIP.rep_keys:
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
            

            # deal with normalization
            if normalize:
                values = values/np.sum(values)
                # TODO: add quantile for density
            else:
                pass
            
            self.value[eCLIP.uID][rep] = values
    
    def window_feature(self, start, end, values):
        """[slice from self.value] 
        start: int, relative position of window
        end: int, relative end position
        sometimes there are window outside of gene length (from point feature)
        """
        if start >=0 and end <= len(values):
            return values[start:end]
        elif start < 0:
            
            valid_values = values[:end]
            pad = [np.nan]*((end-start)-len(valid_values))
            return np.array(pad + list(valid_values))
        elif end > len(values):
            #print('padding to length')
            valid_values = values[start:]
            pad = [np.nan]*((end-start)-len(valid_values))
            return np.array(list(valid_values)+pad)

    def multi_feature_avg(self, intervals, values, align = 'left', max_len = None):
        """[average and align (zero padding) multiple intron/exon
        feature: list of intervals for multiple exons/cds/intron; ex: [(1,100), (2,950)]
        values: 5' to 3' processed value (remove background, normalize) of the whole gene
        align: align feature of different length to the left or right
        return: average feature (np.array)]

        Args:
            intervals ([list]): [list of tuples, each tuple is a interval]
            values ([np.array]): [values of window from]
            align (str, optional): ['left' or 'right', 5 prime=left]. Defaults to 'left'.
            max_len ([int], optional): [description]. Defaults to None.

        Returns:
            [type]: [description]
        """
        all_feature_values = []
        # feature length
        if max_len == None:
            
            
            max_len = max([f[1]-f[0] for f in intervals])
        
        for f in intervals:

            axis = self.to_relative_axis(f)
            feature_values = self.window_feature(start = int(axis[0]),end = int(axis[1]), values = values)
            all_feature_values.append(trim_or_pad(feature_values, max_len, pad_value = np.nan, align = align)) # pad to the same length then append
        
        
        
        feature_average  = np.nanmean(np.stack(all_feature_values), axis = 0)
        return feature_average
        
    def get_average_feature(self, eCLIP, background_method = 'subtract', normalize = True, truncate = True):
        """[fetching values for each feature, also averaging them if needed 
        write values in self.truncate or self.densities, depending on truncate]

        Args:
            eCLIP ([type]): [description]
            background_method (str, optional): [description]. Defaults to 'subtract'.
            normalize (bool, optional): [description]. Defaults to True.
            truncate (bool, optional): [description]. Defaults to True.
        """
        
        # normalize and background
        self.get_value(eCLIP, background_method=background_method, normalize=normalize, truncate=truncate)
        
        # specify saving places
        if truncate:
            to_save = self.truncations  
        else:
            to_save = self.densities

        # initialize empty dictionary to store result
        if eCLIP.uID not in to_save.keys():
            to_save[eCLIP.uID] = {}
        
        for rep in eCLIP.rep_keys:
            if rep not in to_save[eCLIP.uID].keys():
                to_save[eCLIP.uID][rep] = defaultdict(nan_array) # TODO: make good use of defaultdict at the beginning
            values = self.value[eCLIP.uID][rep]
            
            # getting features:
            for feature, fname in zip(self.features.values(), self.features.keys()):
            
                if len(feature) == 0: # no such feature
                    # np.empty does not always yield nan
                    
                    to_save[eCLIP.uID][rep][fname] = nan_array()
                
                elif len(feature) == 1: ### The feature only appears once
                    feature = list(feature)[0]
                    if type(feature) == int: # single point feature
                        feature = (feature-settings.feat_len[fname], feature + settings.feat_len[fname]) # extend a window for it 
                    axis = self.to_relative_axis(feature)
                
                    minus1 = self.window_feature(start = int(axis[0]),end = int(axis[1]), values = values)
                    to_save[eCLIP.uID][rep][fname] = minus1
                
                else: ### The feature only appears multiple time
                    if type(list(feature)[0]) == int:
                        # convert point feature to length
                        feature = [(f-settings.feat_len[fname], f + settings.feat_len[fname])
                                    for f in list(feature)]
                        left_mean = self.multi_feature_avg(feature, values, align = 'left') 
                        to_save[eCLIP.uID][rep][fname] = left_mean

                    else:
                        left_mean = self.multi_feature_avg(feature, values, align = 'left') 
                        right_mean = self.multi_feature_avg(feature, values, align = 'right')
                

                        to_save[eCLIP.uID][rep][fname] = (left_mean, right_mean)
        
 
class Protein_coding_gene(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
        self.featnames+= ['five_prime_UTR','first_CDS', 'CDS', 'last_CDS', 'three_prime_UTR']
    def build_full_cds(self, eCLIP, truncate = True):
        if truncate:
            to_save = self.truncations
        else:
            to_save = self.densities
        
        if eCLIP.uID not in to_save.keys():
            to_save[eCLIP.uID] = {} #TODO: carefully make defaultdict(dict)
        for rep in eCLIP.rep_keys:
            if rep not in to_save[eCLIP.uID].keys():
                to_save[eCLIP.uID][rep] = defaultdict(nan_array)
            to_save[eCLIP.uID][rep]['full_CDS'] = self.concat_multi_feature(self.value[eCLIP.uID][rep], feature = 'CDS')
        
        self.featnames.append('full_CDS') # TODO: make automatic update

        

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
        self.featnames += ['BoxACA', 'BoxC', 'BoxD', 'BoxH', 'Duplex'] 
        
        
class rRNA(Metagene):
    # TODO get from repetitive pipeline
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
class tRNA(Metagene):
    def __init__(self, esnt, chro, start, end, strand, transcript_type, features = None):
        super().__init__(esnt, chro, start, end, strand, transcript_type, features)
        self.featnames += ['anticodon', 't_loop', 'd_loop', 'v_loop', 'anticodon_loop']

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



    
def Build_many_metagene(id_list=None, sample_no = 200, transcript_type = 'protein_coding'):
    """[Create Metagene object for regions for key transcript 
    hg19 does not have 'three_prime_UTR and five_prime_UTR which is annoying]

    Args:
        id_list ([list], optional): [List of transcript/gene IDs]. Defaults to None.
        sample_no (int, optional): [number of genes/transcripts to calculate]. Defaults to 200.
        transcript_type (str, optional): [set to None to build all transcript regardless]. Defaults to 'protein_coding'.
    
    Returns:
        [dict]: [dictionary of ID -> Metagene object]
    """
   # load pre-parsed gencode coords and branchpoint
    # TODO implement None for transcript_type
    # TODO implement canonical transcript filtering
    if transcript_type is None:
        coord_dicts = pickle.load(open(settings.gencode, 'rb'))
    elif transcript_type == 'miRNA':
        coord_dicts = pickle.load(open(settings.mir, 'rb'))
    elif transcript_type == 'snoRNA':
        coord_dicts = pickle.load(open(settings.sno, 'rb'))
    else:
        print('Using:',settings.gencode)
        coord_dicts = pickle.load(open(settings.gencode, 'rb'))
     # Build empty metagene object
    all_metagene = {}
    if id_list==None:
        # if not specify id, use all transcripts
        ids = list(coord_dicts.keys())
        
    elif len(id_list)< sample_no:
        ids = id_list
        
    else: 
        ids = id_list[:sample_no]
    
    for i in ids:
        d = coord_dicts[i]
        if transcript_type is None:
            all_metagene[i] = Metagene.from_dict(d)
        elif d['type'] == transcript_type:
            d['features'] = defaultdict(set, d['features']) # convert to default dict in case feature doesn't exist
            if transcript_type == 'protein_coding':
                all_metagene[i] = Protein_coding_gene.from_dict(d)# TODO: miR, snoR will need another dictionary
            elif transcript_type == 'miRNA':
                all_metagene[i] = miRNA.from_dict(d)
            elif transcript_type == 'tRNA':
                all_metagene[i] = tRNA.from_dict(d)
            elif transcript_type == 'snoRNA':
                all_metagene[i] = snoRNA.from_dict(d)
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



######################################### Relative entropy ##########################################


def density_array_entropy(f, fn):
    '''
    relative entropy for each position, f = positive binding density; fn = negative binding examples's density
    '''
    return np.array([entropy(pk = f[1:, pos], qk = fn[1:, pos]) for pos in range(f.shape[1])])
def nan_array():
            minus1 = np.empty(1)    
            minus1[:] = np.nan
            return minus1

