import pandas as pd
from RBPamp.seed import PSAMSetBuilder
import sys





def get_cluster_set(df, n=100, max_len = 30):
    sorted_cluster = df.sort_values(by = 'pseudoscore', ascending = False)
    kmer_set = []
    i = 0
    while len(kmer_set) < n:
        if len(sorted_cluster.iloc[i]['seq'])< max_len:
            kmer_set.append((sorted_cluster.iloc[i]['seq'], sorted_cluster.iloc[i]['pseudoscore']))
        i += 1
            
        
    return kmer_set



def write_psam_to_file(psams, name):
     for i,p in enumerate(psams):
        
        # save affinity matrix
        df = pd.DataFrame(p.matrix, columns = ['A', 'C', 'G', 'U'])
        df.to_pickle('/home/hsher/encore_region_motif/RBPamp_clus/{}_{}_1000.pickle'.format(name, i))
        
        p.save_logo(fname='/home/hsher/encore_region_motif/RBPamp_clus/{}_{}_1000.svg'.format(name, i))

    

    
if __name__ == '__main__':
    filename = sys.argv[1]
    name = filename.split('/')[-1].replace('.pickle', '')
    df = pd.read_pickle(filename)
    print('loaded {}'.format(name))
    kmer_set = get_cluster_set(df, n = 1000)
    print('got kmer set run RBPamp')
    ps=PSAMSetBuilder(kmer_set)
    psams = ps.make_PSAM_set()
    write_psam_to_file(psams, name)
    print('finish psam, got {} psams'.format(len(psams)))
    
