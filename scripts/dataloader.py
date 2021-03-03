# module that helps to locate all ENCODE eCLIP files fast and easy
# TSCC-specific! Can't use on other computer simply because of the path

import pandas as pd
import os

################## ENCODE 3 ####################
encode_data = pd.read_pickle('~/projects/eclip_encode_id.pickle')
eclip_bam = pd.read_csv('/home/hsher/projects/RBP_annot/ENCODE_FINAL_ANNOTATIONS.uidsonly.txt.manifesthg38.txt', sep = '\t', header= 0)
encode_data = pd.merge(eclip_bam[['uID', 'RBP', 'Cell line']], encode_data, left_on = ['RBP', 'Cell line'], right_on = ['RBP', 'cell_line'])
cols = [e for e in encode_data.columns if 'bam' in e or 'minus' in e or 'plus' in e]
bam_basedir = '/home/hsher/seqdata/eclip_raw/'

for c in cols:
    encode_data[c] = bam_basedir+encode_data[c]
encode_data['idr'] = '/home/hsher/seqdata/eclip_bed/sorted/'+encode_data['uID']+'.01v02.IDR.out.0102merged.bed.blacklist_removed.bed.narrowPeak.bed'
for row in encode_data.index:
    uid = encode_data.loc[row, 'uID']
    encode_data['bed_0'] = '/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20180205/hg38/' + '{0}_0{1}.basedon_{0}_0{1}.peaks.l2inputnormnew.bed.compressed.bed.blacklist_removed.bed'.format(uid, 1)
    encode_data['bed_1'] = '/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20180205/hg38/' + '{0}_0{1}.basedon_{0}_0{1}.peaks.l2inputnormnew.bed.compressed.bed.blacklist_removed.bed'.format(uid, 2)
encode_data['uid'] = encode_data['uID']

################## ENCODE 4 ####################
encode4_data = pd.read_pickle('/home/hsher/projects/ClipNet/ENCODE_stats/ready1214.pickle')
encode4_data['Cell line'] = encode4_data['Cell Line']
master_df = pd.concat([encode_data[['uid','RBP', 'Cell line']], encode4_data[['uid', 'RBP', 'Cell line']]], axis = 0)