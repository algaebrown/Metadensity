
from metadensity.metadensity import *
from metadensity.plotd import *
from metadensity.pos_enrich import *
import pandas as pd
import matplotlib.pyplot as plt
import deepdish as dd
import os

# load encode data
# load IDs
from dataloader import *

outdir='/home/hsher/encore_pos/'
downsampling_path='/home/hsher/subsample_encode/'

def main(uid, sample_no=800, sigma=5, n_largest_to_remove=80):
    ''' run positional enrichment analysis for uid '''
    print('UID:{}'.format(uid), type(uid))
    
    # build eCLIP
    try:
        s = encode_data.loc[encode_data['uid']==str(uid)].iloc[0]
        single_end = False
    except:
        s = encode4_data.loc[encode4_data['uid']==str(uid)].iloc[0]
        single_end = True
    print(single_end, 'Single end')
    
    
    # check if downsampling bam file exists
    bam1=os.path.basename(s['bam_0']).replace('.bam', '.twomillion.bam')
    if os.path.isfile(os.path.join(downsampling_path, bam1)):
        s['bam_0']=os.path.join(downsampling_path, bam1)
        print('using downsampled reads')
    else:
        print('less than two million reads')
    bam2=os.path.basename(s['bam_1']).replace('.bam', '.twomillion.bam')
    if os.path.isfile(os.path.join(downsampling_path, bam2)):
        s['bam_1']=os.path.join(downsampling_path, bam2)

    e = eCLIP.from_series(s, single_end = single_end)

    # build metagene from biogps
    cds_metagenes =  highly_exp_biogps(cell_line = s['Cell line'].upper(),sample_no = sample_no)

    # construct distribution
    m_null, m_ip = construct_distribution(e, cds_metagenes)
    print('Done with {}'.format(e.name))

    print('Calucating enrichment')

    # save raw data
    m_null.save_deepdish(os.path.join(outdir, '{}.null.h5'.format(e.uID)))
    m_ip.save_deepdish(os.path.join(outdir, '{}.ip.h5'.format(e.uID)))
    
    # Wilcoxon
    wx, pval = Wilcox_enrich(m_null, m_ip, n_largest_to_remove = n_largest_to_remove, sigma = sigma) ### need to establish a best practice for it
    dd.io.save(os.path.join(outdir, '{}_wx.pval.h5'.format(e.uID)), pval)
    dd.io.save(os.path.join(outdir, '{}_wx.h5'.format(e.uID)), wx)
    
    # KS test
    ks, pval = KS_enrich(m_null, m_ip, n_largest_to_remove = n_largest_to_remove, sigma = sigma) ### need to establish a best practice for it
    dd.io.save(os.path.join(outdir, '{}_ks.pval.h5'.format(e.uID)), pval)
    dd.io.save(os.path.join(outdir, '{}_ks.h5'.format(e.uID)), ks)
    
if __name__=='__main__':
    uid = sys.argv[1]
    print(uid)
    main(uid)

