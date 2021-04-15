
from metadensity.metadensity import *
from metadensity.plotd import *
from metadensity.pos_enrich import *
import pandas as pd
import matplotlib.pyplot as plt
import deepdish as dd

# load encode data
# load IDs
rbp_file = pd.read_pickle('/home/hsher/projects/ClipNet/ENCODE_stats/rbp_df.pickle')


outdir='/home/hsher/encore_pos/'
def main(uid):
    ''' run positional enrichment analysis for uid '''
    print('UID:{}'.format(uid))
    # build eCLIP
    row = rbp_file.loc[rbp_file['uid']==uid].iloc[0]
    e = eCLIP.from_series(row, single_end = True)

    # build metagene from biogps
    cds_metagenes =  highly_exp_biogps(cell_line = row['Cell Line'].upper(),sample_no = 100)

    # construct distribution
    m_null, m_ip = construct_distribution(e, cds_metagenes)
    print('Done with {}'.format(e.name))

    print('Calucating enrichment')
    # run AUC
    #auc = AUC_enrich(m_null, m_ip)
    #ks, pval = AUC_enrich(m_null, m_ip)

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

