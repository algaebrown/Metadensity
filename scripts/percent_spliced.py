import pysam
########### JUNCTION READ RELATED ##############
def is_junction_read(cigartuple, intron_len = 500):
    Ns = [t for t in cigartuple if t[0]==3]
    if len(Ns)>0 and Ns[0][1]>intron_len:
        return True
    else:
        return False

def percent_junction_read(reads):
    ''' return % of reads are junctional'''
    if len(reads)==0:
        return 0
    else:
        perc_junc = sum([is_junction_read(r.cigartuples) for r in reads])/len(reads)
        
    
        return perc_junc

import sys
sys.path.append('/home/hsher/Metadensity/scripts')
import pandas as pd
import matplotlib.pyplot as plt
from dataloader import *
from metadensity.truncation import *
from scipy.stats import ttest_rel
def perc_junctional_all(bams, single_end = False, n_read = 0, n_transcript = 1000):
    transcripts = BedTool('/home/hsher/gencode_coords/gencode.v33.transcript.gff3')
    data = []
    for t in transcripts[:n_transcript]:
        row = [t.attrs['ID']]
        for bam in bams:

            reads = fetch_reads(bam, chrom = t.chrom,
                                start = t.start,
                                end = t.end,
                                strand = t.strand,
                                single_end = single_end)

            if len(reads)<=n_read:
                row.append(None)
            else:

                p=percent_junction_read(reads)
                row.append(p)
        data.append(row)
    if len(bams)==3:
        df = pd.DataFrame(data, columns = ['tid','bam1', 'bam2', 'bamin']).dropna()
    else:
        
        df = pd.DataFrame(data, columns = ['tid','bam1', 'bam2', 'bamin1', 'bamin2']).dropna()
    return df
        
def test_junctional_reads(uid, n_read = 10, n_transcript = 10000):
    try:
        bams= return_fobj3(uid)
        single_end = False
    except:
        bams = return_fobj4(uid)
        single_end = True
    
    data = perc_junctional_all(bams, single_end = single_end, n_read = n_read, n_transcript = n_transcript)
    data.to_csv(f'/home/hsher/densities/{uid}.perc_spliced.csv')
    
    # testing: more junctional reads in IP
    if len(bams) == 3:
        tstat1, gr_pv1 = ttest_rel(data['bam1'], data['bamin'], alternative = 'greater')
        tstat2, gr_pv2 = ttest_rel(data['bam2'], data['bamin'], alternative = 'greater')
    else:
        tstat1, gr_pv1 = ttest_rel(data['bam1'], data['bamin1'], alternative = 'greater')
        tstat2, gr_pv2 = ttest_rel(data['bam2'], data['bamin2'], alternative = 'greater')
        
    if len(bams) == 3:
        tstat1, ls_pv1 = ttest_rel(data['bam1'], data['bamin'], alternative = 'less')
        tstat2, ls_pv2 = ttest_rel(data['bam2'], data['bamin'], alternative = 'less')
    else:
        
        tstat1, ls_pv1 = ttest_rel(data['bam1'], data['bamin1'], alternative = 'less')
        tstat2, ls_pv2 = ttest_rel(data['bam2'], data['bamin2'], alternative = 'less')
        
    return gr_pv1, gr_pv2, ls_pv1, ls_pv2
    
    
    
all_data = []
for uid in master_df['uid']:
    try:
        gr_pv1, gr_pv2, ls_pv1, ls_pv2 = test_junctional_reads(uid)
        all_data.append([uid, gr_pv1, gr_pv2, ls_pv1, ls_pv2])
    except Exception as e:
        print(uid, e)

all_data = pd.DataFrame(all_data, columns = ['uid', 'greater_pval_rep1', 'greater_pval_rep2',
                                            'lesser_pval_rep1','lesser_pval_rep2'])
all_data['RBP'] = all_data['uid'].map(master_df.set_index('uid')['RBP'])
from statsmodels.stats.multitest import fdrcorrection
for col in [col for col in all_data.columns if 'pval' in col]:
    _, all_data[col.replace('pval', 'fdr')] = fdrcorrection(all_data[col])

all_data.to_csv('~/projects/Metadensity/notebooks/perc_splice.csv')