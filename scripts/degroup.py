from metadensity.metadensity import *
from metadensity.plotd import *
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import sys
from dataloader import *
from metadensity.pos_enrich import *

outdir='/home/hsher/de_metaden/'
# parse DE filenames
de_path = '/home/hsher/deseq_gccor/normalized/'
all_de_files = os.listdir(de_path)
all_de_files.remove('result_URLs_HepG2.txt')
all_de_files.remove('result_URLs_K562.txt')
all_de_files_rbp = [f.split('-')[0] for f in all_de_files if '-' in f]
all_de_files_cell = [f.split('-')[2].split('_')[0] for f in all_de_files if '-' in f]
de_df = pd.DataFrame([all_de_files, all_de_files_rbp, all_de_files_cell]).T
de_df.columns = ['fname', 'rbp', 'cell line']

def DEseq_result(fname):
    ''' return transcript ID that is up and down regulated'''
    # read DEseq results
    de = pd.read_csv(de_path+ fname, sep = '\t', header = 0, index_col = 0)
    
    # select those significant
    de = de[de['padj']< 0.05]
    de.reset_index(inplace = True)
    
    # select fold change high
    up = de.loc[de['log2FoldChange']> 1]
    down = de.loc[de['log2FoldChange']< -1]
    
    # remove version
    up = [i.split('.')[0] for i in up['index'].tolist()]
    down = [i.split('.')[0] for i in down['index'].tolist()]
    return up, down
    
def find_direct_bind(id_list, eCLIP, filter_direct = True, enough_transcript = None):
    ''' find DE genes that has at least 1 idr peak'''
    subset = transcript.filter(lambda x: x.attrs['gene_id'].split('.')[0] in id_list and x.attrs['ID'] in enough_transcript).saveas()
    if filter_direct:
        direct_bind_transcript = subset.intersect(eCLIP.idr, s = True, u = True).saveas()
    else:
        return subset
    
    return direct_bind_transcript

def get_de_metagroup(cell_line, rbp, filter_direct = True, thres = 100):
    # start eCLIP object for you
    try:
        s = encode_data.loc[(encode_data['RBP'] == rbp)&(encode_data['Cell line'] == cell_line)].iloc[0]
        e = eCLIP.from_series(s)
    except:
        s = encode4_data.loc[(encode4_data['RBP'] == rbp)&(encode4_data['Cell line'] == cell_line)].iloc[0]
        e = eCLIP.from_series(s, single_end = True)
        print('using encode4')

        
    transcripts_with_reads = e.enough_transcripts(thres = thres)
    print('NUMBER OF TRANSCRIPT WITH ENOUGH:', len(transcripts_with_reads))
    
    # fetch DE genes
    fname = de_df.loc[(de_df['rbp'] == rbp) & (de_df['cell line'] == cell_line), 'fname'].values[0]
    up, down = DEseq_result(fname)
    
    up_transcripts = find_direct_bind(up, e, filter_direct = filter_direct, enough_transcript = transcripts_with_reads)
    down_transcripts = find_direct_bind(down, e, filter_direct = filter_direct, enough_transcript = transcripts_with_reads)
    
    backgrounds = set(transcripts_with_reads) - set([t.attrs['ID'] for t in up_transcripts])- set([t.attrs['ID'] for t in down_transcripts])
    
    print('# up reg: {} # down reg: {} background: {}'.format(len(up_transcripts), len(down_transcripts), len(backgrounds)))

    meta_back = Metadensity(e, name=rbp+'_nochange', transcript_ids= list(backgrounds), background_method = 'relative information', normalize = False, sample_no = 500)
    meta_back.get_density_array()
    if len(up_transcripts) > 0:
        meta_up = Metadensity(e, name=rbp+'_upreg', transcripts = up_transcripts, background_method = 'relative information', normalize = False, sample_no = len(up_transcripts))
        meta_up.get_density_array()
    else:
        meta_up=None
    
    if len(down_transcripts) > 0:
        meta_down = Metadensity(e, name=rbp+'_downreg', transcripts = down_transcripts, background_method = 'relative information', normalize = False, sample_no = len(down_transcripts))
        meta_down.get_density_array()
    else:
        meta_down=None
    
    # then we create metadensity for each of them
    
    return [meta_up, meta_down, meta_back]

def main(cell_line, rbp, filter_direct = False, thres = 100, sigma=5):
    ''' get metadensity for each group of protein, do stats '''
    metas = get_de_metagroup(cell_line, rbp, filter_direct = filter_direct, thres = thres) # up, down, no change group

    

    for m,group in zip(metas, ['up', 'down', 'nc']):
        if m is not None:
            # sometimes has no up/down transcripts
            m.save_deepdish(os.path.join(outdir, f'{rbp}_{cell_line}.{group}.denarray.h5'))
    
    if metas[0] is not None:
        ks_up,p_up=KS_enrich(metas[-1], metas[0], bidir=True, sigma = sigma)
        dd.io.save(os.path.join(outdir, f'{rbp}_{cell_line}.up.ks.h5'), ks_up)
        dd.io.save(os.path.join(outdir, f'{rbp}_{cell_line}.up.ks.pval.h5'), p_up)

    if metas[1] is not None:
        ks_down,p_down=KS_enrich(metas[-1], metas[1], bidir=True, sigma = sigma)
        dd.io.save(os.path.join(outdir, f'{rbp}_{cell_line}.down.ks.h5'), ks_down)
        dd.io.save(os.path.join(outdir, f'{rbp}_{cell_line}.down.ks.pval.h5'), p_down)

if __name__=='__main__':
    rbp=sys.argv[1]
    cell_line=sys.argv[2]

    main(cell_line, rbp)