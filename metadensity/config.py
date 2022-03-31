# set up genome coordinates
import os
from collections import defaultdict
import configparser

# change function
class MetadensityConfig:

    def __init__(self):
        print('please set the right config according to genome coordinate')
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