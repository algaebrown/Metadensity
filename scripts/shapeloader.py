import os
import numpy as np
import pyBigWig
encode_shape_dir = '/home/hsher/icshape_data'
zhang_shape_dir = '/home/hsher/icshape_data_zhang'
shape_data_dir = {'ENCODE':{'K562':os.path.join(encode_shape_dir, 'ENCFF161UHZ.tsv'),
                            'HepG2': os.path.join(encode_shape_dir, 'ENCFF706GXJ.tsv')
                            },
                    'Zhang':{'K562':os.path.join(zhang_shape_dir, 'GSM4333260_K562.out.txt'),
                            'HepG2':os.path.join(zhang_shape_dir, 'GSM4333262_HepG2.out.txt')
                            }
                }

class SHAPE_data:
    def __init__(self, cell_line):
        self.plus = pyBigWig.open('/home/hsher/icshape_data_zhang/icSHAPE/{}-plus.bw'.format(cell_line))
        self.minus = pyBigWig.open('/home/hsher/icshape_data_zhang/icSHAPE/{}-minus.bw'.format(cell_line))
        self.cell_line = cell_line
    def fetch(self, chrom = None, start= None, end=None, strand= None, interval = None):
        ''' return icSHAPE reacitivity for a bedtool interval or chrom, start, end, strand'''
        if interval:
            start = interval.start
            end = interval.end
            strand = interval.strand
            chrom = interval.chrom
        if strand == '-':
            icshape_data = self.minus
        else:
            icshape_data = self.plus
        values = icshape_data.values(chrom, start, end)
        if strand == '-':
            values = values[::-1]
        return values

class fSHAPE_data(SHAPE_data):
    def __init__(self, cell_line, treatment, rep):
        self.plus = pyBigWig.open(f'/home/hsher/fshape_raw/CITS/{cell_line}_{treatment}{rep}.plus.bw')
        self.minus = pyBigWig.open(f'/home/hsher/fshape_raw/CITS/{cell_line}_{treatment}{rep}.minus.bw')
        self.cell_line = cell_line
        self.treatment = treatment
        self.rep = rep


def read_icshape(fname):
    ''' read icSHAPE data from tsv file 
    return dict, keys: transcript_id or gene_id; values: 5' to 3' icSHAPE reactivity in list'''
    data = {}
    with open(fname) as f:
        for line in f:
            
            values = line[:-1].split('\t') # remove \n
            gene_id = values[0]
            gene_len = values[1]
            coverage = values[2]
            reactivity = [float(v) if v!= 'NULL' else np.nan for v in values[3:]]
            data[gene_id] = reactivity
    return data