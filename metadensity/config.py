# set up genome coordinates
import os
from collections import defaultdict
import configparser

# change function
class MetadensityConfig:

    def __init__(self):
        # original file
        self.root_dir = '/home/hsher/gencode_coords/'
        self.genome_coord = 'hg38'

        self.file_dict = {
            'hg38': {
                'genome':'GRCh38.p13.genome.fa',
                'gencode':'gencode.v33.annotation.gff3',
                'transcript':'gencode.v33.transcript.gff3', # transcripts, canonical 
                'gencode_feature': 'gencode.v33.combine.sorted.gff3',
                'branchpoint': 'branchpoint_hg38.bed', # experimental branchpoint
                'branchpoint_pred': 'gencode_v26_branchpoints_pred.csv', # predicted branchpoint
                'polyA': '/home/hsher/gencode_coords/polyA.atlas.clusters.2.0.GRCh38.96.bed',
                'miR': 'miR/hsa_hg38.gff3',
                'pri-mir': '', #TODO lifeover
                'snoRNA': 'snoDB_hg38.tsv', # TODO call boxes
                'tRNA': '', # TODO call loops
                'lncRNA':'lncipedia_5_2_hc_hg38.gff'
            },
            'hg19':{
                'genome':'',
                'gencode': 'gencode.v19.annotation.gff3',
                'transcript': 'gencode.v19.transcript.gff3',
                'gencode_feature': 'gencode.v19.combine.sorted.gff3',
                'branchpoint': 'branchpoint_hg19.bed',
                'branchpoint_pred': 'gencode_v19_branchpoints_pred.csv',
                'polyA': '', #TODO liftover
                'miR': 'miR/hsa_hg19.gff3',
                'pri-mir':'human_primiR_hg19.bed',
                'snoRNA': '', # TODO liftover
                'tRNA':'',
                'lncRNA':'lncipedia_5_2_hc_hg19.gff'
            },
            'mm10':{
                'gencode': 'gencode.vM25.annotation.gff3',
                'transcript': 'gencode.vM25.transcript.gff3',
                'gencode_feature': 'gencode.vM25.combine.sorted.gff3',
                'branchpoint': '',
                'branchpoint_pred': '',
                'polyA': '', #TODO liftover
                'miR': '',
                'pri-mir':'',
                'snoRNA': '', # TODO liftover
                'tRNA':'',
                'lncRNA':''
            },
            'ChlSab':{
                'gencode': 'Chlor.gtf',
                'transcript': 'Chlor.genes.gtf',
                'gencode_feature': 'Chlor.combined.sorted.gtf',
                'branchpoint': '',
                'branchpoint_pred': '',
                'polyA': '', #TODO liftover
                'miR': '',
                'pri-mir':'',
                'snoRNA': '', # TODO liftover
                'tRNA':'',
                'lncRNA':''
                },
                
        }

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
        self.refresh_coords()
    def refresh_coords(self):
        self.fasta = os.path.join(self.root_dir, self.file_dict[self.genome_coord]['genome'])
        self.gencode_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['gencode'])
        self.transcript_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['transcript'])
        self.gencode_feature_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['gencode_feature'])
        
        self.branchpoint_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['branchpoint'])
        self.branchpoint_pred_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['branchpoint_pred'])
        self.polya_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['polyA'])
        self.miR_feature_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['miR'])
        self.pri_miR_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['pri-mir'])
        self.snoRNA_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['snoRNA'])
        self.tRNA_fname = os.path.join(self.root_dir,self.file_dict[self.genome_coord]['tRNA'])
        self.lncRNA_fname= os.path.join(self.root_dir,self.file_dict[self.genome_coord]['lncRNA'])
        
        # processed coords
        self.datadir = os.path.join(os.path.dirname(__file__), 'data', self.genome_coord)
        self.gencode = os.path.join(self.datadir, 'gencode')
        self.mir = os.path.join(self.datadir, 'miRNA')
        self.sno = os.path.join(self.datadir, 'snoRNA')
        self.lnc = os.path.join(self.datadir, 'lnc')
        self.trna = os.path.join(self.datadir, 'tRNA')

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