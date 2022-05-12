# set up genome coordinates
import os
from collections import defaultdict
import configparser

# change function
class MetadensityConfig:

    def __init__(self):
        print('please set the right config according to genome coordinate')
        self.from_config_file('/opt/Metadensity/config/hg38.ini')
        print(f'Using HG38 by default')

        self.feat_len = {'five_prime_UTR':100,  ####### generic DNA
                        'three_prime_UTR':150, 
                        'intron':500,
                        'exon':100,
                        'first_exon':100,
                        'last_exon': 100,
                        'CDS': 100, ######## protein coding
                        'first_CDS': 100,
                        'last_CDS':100,
                        '3p_duplex': 20, ############ miR
                        '5p_duplex': 20,
                        'mature': 20,
                        'hairpin':10,
                        'pri_mir': 500, 
                        'start_codon':3,
                        'stop_codon':3,
                        'branchpoint':50,
                        'branchpoint_pred':50,
                        'full_CDS':200}
        self.point_feature_len=50
        self.feat_len = defaultdict(lambda: self.point_feature_len, self.feat_len)
        
        
        self.ax_width_dict = {'UTR':1,
                'exon':2,
                'intron':3,
                'CDS':2,
                'duplex':2,
                'mature':2,
                'hairpin':2,
                'branchpoint':2}
        self.ax_width_dict = defaultdict(lambda: 2, self.ax_width_dict)

    def from_config_file(self, filename):
        cps = configparser.ConfigParser()
        cps.read(filename)

        self.fasta = cps['FILES']['GENOME_FA']
        self.gencode_fname = cps['FILES']['GENCODE']
        self.transcript_fname = cps['FILES']['TRANSCRIPT']
        self.gencode_feature_fname = cps['FILES']['FEATURE']

        print(f'Using {self.fasta}')
        
        self.branchpoint_fname = cps['FILES']['BRANCHPOINT']
        self.branchpoint_pred_fname = cps['FILES']['BRANCHPOINT_PRED']
        self.polya_fname = cps['FILES']['POLYA']
        self.miR_feature_fname =  cps['FILES']['MIRNA']
        self.snoRNA_fname = cps['FILES']['SNORNA']
        self.lncRNA_fname= cps['FILES']['LNCRNA']
        
        # processed coords
        self.datadir = cps['FILES']['DATADIR']
        self.gencode = os.path.join(self.datadir, 'gencode')
        self.mir = os.path.join(self.datadir, 'miRNA')
        self.sno = os.path.join(self.datadir, 'snoRNA')
        self.lnc = os.path.join(self.datadir, 'lnc')
        self.trna = os.path.join(self.datadir, 'tRNA')

        

        


settings = MetadensityConfig()
