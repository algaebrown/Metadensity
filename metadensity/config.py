# set up genome coordinates


# change function
class MetadensityConfig:

    def __init__(self):
        self.root_dir = '/home/hsher/gencode_coords/'
        self.transcript_fname = 'gencode.v33.transcript.gff3'
        self.gencode_feature_fname = 'gencode.v33.combine.sorted.gff3'
        self.miR_feature_fname = 'miR/hsa_hg38.gff3'

        self.feat_len = {'five_utr':100,  ####### generic DNA
                        'three_utr':150, 
                        'intron':1500,
                        'exon':150,
                        'first_exon':150,
                        'last_exon': 150,
                        'CDS': 150, ######## protein coding
                        'first_CDS': 150,
                        'last_CDS':150,
                        '3p_duplex': 20, ############ miR
                        '5p_duplex': 20,
                        'mature': 20,
                        'hairpin':10,
                        'pri_mir': 500}
        
        self.ax_width_dict = {'utr':1,
                'exon':2,
                'intron':3,
                'CDS':2,
                'duplex':2,
                'mature':2,
                'hairpin':2}


settings = MetadensityConfig()