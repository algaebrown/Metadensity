# module that helps to locate all ENCODE eCLIP files fast and easy
# TSCC-specific! Can't use on other computer simply because of the path

import pandas as pd
import os
import pysam
from dataloader import *

################## ENCODE 3 ####################
encode_data = pd.read_pickle('/home/hsher/projects/ClipNet/ENCODE_stats/eclip_encode_id.pickle')
################## ENCODE 4 ####################
encode4_data =pd.read_pickle('/home/hsher/projects/ClipNet/ENCODE_stats/rbp_df.pickle')

master_df = pd.concat([encode_data, encode4_data], axis = 0)


############ load bam ###################
def return_fobj3(uid, encode_data = encode_data):
    '''return BedTool, Pysam objects'''
    
    
    row = encode_data.loc[encode_data['uid']==uid]
    if row.shape[0] == 0:
        print('No matching data in ENCODE 3')
    else:
        bam1 = row['bam_0'].values[0]
        bam2 = row['bam_1'].values[0]
        bam_in = row['bam_control'].values[0]
    
    
    bam1_fobj = pysam.Samfile(bam1, 'rb')
    bam2_fobj = pysam.Samfile(bam2, 'rb')
    bam_input_fobj = pysam.Samfile(bam_in, 'rb')
    
    return bam1_fobj, bam2_fobj, bam_input_fobj
def return_fobj4(uid, rbp_file = encode4_data):
    '''return BedTool, Pysam objects'''
    
    
    row = rbp_file.loc[rbp_file['uid']==uid]
    if row.shape[0] == 0:
        print('No matching data in ENCODE 4')
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