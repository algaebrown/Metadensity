import pandas as pd
from RBPamp.seed import PSAMSetBuilder

z_df_filt_95 = pd.read_pickle('~/dbclus_motif/encore_zscore_filt95.pickle')

def get_kmer_set(rbp):
    ''' from kmer zscore df, return kmer_set for PSAM_builder'''
    s = z_df_filt_95.loc[rbp].sort_values(ascending = False)
    kmer_set = []
    for kmer in s.index:
        kmer_set.append((kmer, s[kmer]))
    return kmer_set


def run_rbpamp(index):
    kmer_set = get_kmer_set(index)
    ps=PSAMSetBuilder(kmer_set[:100])
    psams = ps.make_PSAM_set()
    
    for i,p in enumerate(psams):
        
        # save affinity matrix
        df = pd.DataFrame(p.matrix, columns = ['A', 'C', 'G', 'U'])
        df.to_pickle('/home/hsher/encore_region_motif/RBPamp_kmer/{}_{}.pickle'.format(index, i))
        
        p.save_logo(fname='/home/hsher/encore_region_motif/RBPamp_kmer/{}_{}.svg'.format(index, i))

for index in z_df_filt_95.index:
    run_rbpamp(index)
    print('done with {}'.format(index))