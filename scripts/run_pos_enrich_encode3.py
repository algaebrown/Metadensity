
from metadensity.metadensity import *
from metadensity.plotd import *
from metadensity.pos_enrich import *
import pandas as pd
import matplotlib.pyplot as plt
import deepdish as dd

# load encode data
# load IDs
encode_data = pd.read_pickle('~/projects/eclip_encode_id.pickle')
eclip_bam = pd.read_csv('/home/hsher/projects/RBP_annot/ENCODE_FINAL_ANNOTATIONS.uidsonly.txt.manifesthg38.txt', sep = '\t', header= 0)

# join data
encode_data = pd.merge(eclip_bam[['uID', 'RBP', 'Cell line']], encode_data, left_on = ['RBP', 'Cell line'], right_on = ['RBP', 'cell_line'])

cols = [e for e in encode_data.columns if 'bam' in e or 'minus' in e or 'plus' in e]
bam_basedir = '/home/hsher/seqdata/eclip_raw/'
import os
for c in cols:
    encode_data[c] = bam_basedir+encode_data[c]
encode_data['idr'] = '/home/hsher/seqdata/eclip_bed/sorted/'+encode_data['uID']+'.01v02.IDR.out.0102merged.bed.blacklist_removed.bed.narrowPeak.bed'
for row in encode_data.index:
    uid = encode_data.loc[row, 'uID']
    encode_data['bed_0'] = '/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20180205/hg38/' + '{0}_0{1}.basedon_{0}_0{1}.peaks.l2inputnormnew.bed.compressed.bed.blacklist_removed.bed'.format(uid, 1)
    encode_data['bed_1'] = '/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20180205/hg38/' + '{0}_0{1}.basedon_{0}_0{1}.peaks.l2inputnormnew.bed.compressed.bed.blacklist_removed.bed'.format(uid, 2)
encode_data['uid'] = encode_data['uID']


outdir='/home/hsher/encore_pos/'
def main(uid):
    ''' run positional enrichment analysis for uid '''
    print('UID:{}'.format(uid), type(uid))
    # build eCLIP
    s = encode_data.loc[encode_data['uID']==str(uid)].iloc[0]
    print(s)
    e = eCLIP.from_series(s, single_end = False)

    # build metagene from biogps
    cds_metagenes =  highly_exp_biogps(cell_line = s['cell_line'].upper(),sample_no = 100)

    # construct distribution
    m_null, m_ip = construct_distribution(e, cds_metagenes)
    print('Done with {}'.format(e.name))

    print('Calucating enrichment')
    # run AUC
    #ks, pval = AUC_enrich(m_null, m_ip)

    # run reject
    #rej = enrich_by_thres(m_null, m_ip)
    wx, pval = Wilcox_enrich(m_null, m_ip, n_largest_to_remove = 20) ### need to establish a best practice for it

    # run reject
    #rej = enrich_by_thres(m_null, m_ip)

    print('saving results')
    # save result
    #dd.io.save(os.path.join(outdir, '{}_auc.h5'.format(e.uID)), auc)
    #dd.io.save(os.path.join(outdir, '{}_ks.h5'.format(e.uID)), ks)
    dd.io.save(os.path.join(outdir, '{}_pval.h5'.format(e.uID)), pval)
    dd.io.save(os.path.join(outdir, '{}_wx.h5'.format(e.uID)), wx)
    

if __name__=='__main__':
    uid = sys.argv[1]
    print(uid)
    main(uid)

