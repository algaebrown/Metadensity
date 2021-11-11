import pandas as pd
import sys
import numpy as np
from RBPamp.seed import PSAMSetBuilder

def get_kmer_set(z_df, sort_col = 'zscore', ascending = False):
    ''' from kmer zscore df, return kmer_set for PSAM_builder'''
    kmer_series = z_df[sort_col].sort_values(ascending = ascending).dropna()
    kmer_set = []
    for kmer in kmer_series.index:
        kmer_set.append((kmer, kmer_series[kmer]))
    return kmer_set


def run_rbpamp(z_df, out_prefix, sort_col = 'zscore', ascending = False, top_kmer = 100):
    kmer_set = get_kmer_set(z_df, sort_col = sort_col, ascending = ascending)
    ps=PSAMSetBuilder(kmer_set[:top_kmer])

    # generating multiple PSAM
    psams = ps.make_PSAM_set()
    
    for i,p in enumerate(psams):
        
        # save affinity matrix
        df = pd.DataFrame(p.matrix, columns = ['A', 'C', 'G', 'U'])
        df.to_csv('{}.psam.{}.{}.csv'.format(out_prefix, sort_col, i))
        print('==== PSAM {} ====='.format(i))
        print(df)
        print('==================')
        
        # save figure
        p.save_logo(fname='{}.psam.{}.{}.svg'.format(out_prefix, sort_col, i))
    return psams



if __name__ == "__main__":
    # load z-scores
    outfile = sys.argv[1] # zscore, pseudoscore_ks_pval, csln_entropy
    out_prefix = outfile.replace('.csv', '')
    #outfile = os.path.join('/home/hsher/dbclus_motif', prefix+'{}.{}mer.csv'.format(rep, k))

    # 
    z_df = pd.read_csv(outfile, index_col = 0, header = 0)
    z_df['logp'] = -np.log10(z_df['pseudoscore_ks_pval'])
    run_rbpamp(z_df, out_prefix, sort_col = 'zscore', ascending = False)
    run_rbpamp(z_df, out_prefix, sort_col = 'logp', ascending = False)