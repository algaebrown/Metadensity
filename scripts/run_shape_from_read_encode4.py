from metadensity.shape_from_read import *
from pybedtools import BedTool
import pandas as pd
import pysam
import os
import numpy as np
import sys


rbp_file = pd.read_pickle('/home/hsher/projects/ClipNet/ENCODE_stats/rbp_df.pickle')
def return_fobj(uid, rbp_file = rbp_file):
    '''return BedTool, Pysam objects'''
    raw_path = '/home/hsher/seqdata/eclip_raw/'
    
    row = rbp_file.loc[rbp_file['uid']==uid]
    if row.shape[0] == 0:
        print('No matching data')
    else:
        bam1 = row['bam_0'].values[0]
        bam2 = row['bam_1'].values[0]
        bam_in1 = row['bam_control_0'].values[0]
        bam_in2 = row['bam_control_1'].values[0]
    
    
    bam1_fobj = pysam.Samfile(bam1, 'rb')
    bam2_fobj = pysam.Samfile(bam2, 'rb')
    bam_input1_fobj = pysam.Samfile(bam_in1, 'rb')
    bam_input2_fobj = pysam.Samfile(bam_in2, 'rb')
    
    return bam1_fobj, bam2_fobj, bam_input1_fobj, bam_input2_fobj

def main(uid):
    # read bam files
    bam1, bam2, bam_in1, bam_in2 = return_fobj(uid)

    # cell line
    cell_line = rbp_file.loc[rbp_file['uid']==uid, 'Cell Line'].iloc[0]
    #cell_line = encode_data.loc[encode_data['uID']==uid, 'Cell line'].iloc[0]

    # load SHAPE data as is
    print('loading SHAPE for cell line {}'.format(cell_line))
    shape_dir = '/home/hsher/icshape_data'
    if cell_line == 'HepG2':
        
        data = read_icshape(os.path.join(shape_dir, 'ENCFF706GXJ.tsv'))
    else:
        data = read_icshape(os.path.join(shape_dir, 'ENCFF161UHZ.tsv'))
    
    # filter for features
    coords = BedTool('/home/hsher/gencode_coords/gencode.v29.gene.gff3')
    # filter for genes with icSHAPE data
    filtered = coords.filter(lambda x: x.attrs['ID'] in data.keys()).saveas()
    count_overlap = filtered.merge(s = True, c = [1,6,7], o = ['count','distinct','distinct']).saveas()
    is_overlap = count_overlap.filter(lambda x: int(x[3])>1).saveas()
    filtered = filtered.intersect(is_overlap, s = True, v = True)

    # get shape values around read start sites
    w_ip1, w_ip2, w_in1, w_in2 = windows_for_all(bam1, bam2, bam_in1,  data, filtered,  bam_in2=bam_in2, window = 50,  num = len(filtered), single_end = True)
    
    # calculate enrichment
    ps1, ks1 = ks_all_pos(w_ip1, w_in1)
    ps2, ks2 = ks_all_pos(w_ip2, w_in2)

    outdir = '/home/hsher/eclip_read_shape'
    
    np.save(os.path.join(outdir, '{}_ip1'.format(uid)), w_ip1)
    np.save(os.path.join(outdir, '{}_ip2'.format(uid)), w_ip2)
    np.save(os.path.join(outdir, '{}_in1'.format(uid)), w_in1)
    np.save(os.path.join(outdir, '{}_in2'.format(uid)), w_in2)

    print(w_ip1.shape, w_ip2.shape, w_in1.shape, w_in2.shape)
    
    np.save(os.path.join(outdir, '{}_1.pval'.format(uid)), ps1)
    np.save(os.path.join(outdir, '{}_2.pval'.format(uid)), ps2)

    np.save(os.path.join(outdir, '{}_1.ks'.format(uid)), ks1)
    np.save(os.path.join(outdir, '{}_2.ks'.format(uid)), ks2)

if __name__=='__main__':
    uid = sys.argv[1]
    main(uid)