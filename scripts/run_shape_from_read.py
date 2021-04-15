from metadensity.shape_from_read import *
from pybedtools import BedTool
import pandas as pd
import pysam
import os
import numpy as np
import sys

encode_data = pd.read_pickle('~/projects/eclip_encode_id.pickle')
eclip_bam = pd.read_csv('/home/hsher/projects/RBP_annot/ENCODE_FINAL_ANNOTATIONS.uidsonly.txt.manifesthg38.txt', sep = '\t', header= 0)

# join data
encode_data = pd.merge(eclip_bam[['uID', 'RBP', 'Cell line']], encode_data, left_on = ['RBP', 'Cell line'], right_on = ['RBP', 'cell_line'])
def return_fobj(uid, encode_data = encode_data):
    '''return BedTool, Pysam objects'''
    raw_path = '/home/hsher/seqdata/eclip_raw/'
    
    row = encode_data.loc[encode_data['uID']==uid]
    if row.shape[0] == 0:
        print('No matching data')
    else:
        bam1 = row['bam_0'].values[0]
        bam2 = row['bam_1'].values[0]
        bam_in = row['bam_control'].values[0]
    
    
    bam1_fobj = pysam.Samfile(raw_path + bam1, 'rb')
    bam2_fobj = pysam.Samfile(raw_path + bam2, 'rb')
    bam_input_fobj = pysam.Samfile(raw_path + bam_in, 'rb')
    
    return bam1_fobj, bam2_fobj, bam_input_fobj

def main(uid):
    # read bam files
    bam1, bam2, bam_in = return_fobj(uid)

    # cell line
    cell_line = encode_data.loc[encode_data['uID']==uid, 'Cell line'].iloc[0]

    # load SHAPE data as is
    print('loading SHAPE for cell line {}'.format(cell_line))
    shape_dir = '/home/hsher/icshape_data'
    if cell_line == 'HepG2':
        
        data = read_icshape(os.path.join(shape_dir, 'ENCFF706GXJ.tsv'))
    else:
        data = read_icshape(os.path.join(shape_dir, 'ENCFF161UHZ.tsv'))
    
    # filter for features
    coords = BedTool('/home/hsher/gencode_coords/gencode.v29.gene.gff3')
    # filter for genes with icSHAPE data, make sure the regions are not overlapping so that we don't take some data twice.
    filtered = coords.filter(lambda x: x.attrs['ID'] in data.keys()).saveas()
    count_overlap = filtered.merge(s = True, c = [1,6,7], o = ['count','distinct','distinct']).saveas()
    is_overlap = count_overlap.filter(lambda x: int(x[3])>1).saveas()
    filtered = filtered.intersect(is_overlap, s = True, v = True)

    # get shape values around read start sites
    w_ip1, w_ip2, w_in = windows_for_all(bam1, bam2, bam_in, data, filtered, bam_in2 = None, window = 50,  num = 100, single_end = False)
    
    # calculate enrichment
    ps1, ks1 = ks_all_pos(w_ip1, w_in)
    ps2, ks2 = ks_all_pos(w_ip2, w_in)

    outdir = '/home/hsher/eclip_read_g4'
    
    np.save(os.path.join(outdir, '{}_ip1'.format(uid)), w_ip1)
    np.save(os.path.join(outdir, '{}_ip2'.format(uid)), w_ip2)
    np.save(os.path.join(outdir, '{}_in'.format(uid)), w_in)
    
    np.save(os.path.join(outdir, '{}_1.pval'.format(uid)), ps1)
    np.save(os.path.join(outdir, '{}_2.pval'.format(uid)), ps2)

    np.save(os.path.join(outdir, '{}_1.ks'.format(uid)), ks1)
    np.save(os.path.join(outdir, '{}_2.ks'.format(uid)), ks2)

if __name__=='__main__':
    uid = sys.argv[1]
    main(uid)