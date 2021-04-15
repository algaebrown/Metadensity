import os
import numpy as np
encode_shape_dir = '/home/hsher/icshape_data'
zhang_shape_dir = '/home/hsher/icshape_data_zhang'
shape_data_dir = {'ENCODE':{'K562':os.path.join(encode_shape_dir, 'ENCFF161UHZ.tsv'),
                            'HepG2': os.path.join(encode_shape_dir, 'ENCFF706GXJ.tsv')
                            },
                    'Zhang':{'K562':os.path.join(zhang_shape_dir, 'GSM4333260_K562.out.txt'),
                            'HepG2':os.path.join(zhang_shape_dir, 'GSM4333262_HepG2.out.txt')
                            }
                }


def read_icshape(fname):
    ''' read icSHAPE data from file 
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