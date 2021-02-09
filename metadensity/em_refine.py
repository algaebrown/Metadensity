import sys
sys.path.append('/home/hsher/projects/Metadensity/metadensity')
import os
from seqem import *
from sequence import simulate_kmer_background, kmer_zscore
import seaborn as sns
import logging
from optparse import OptionParser

## input: control_df, df
## conda environment: rbpamp


def main(df, ctrl_df, k = 7, niters = 10, pseudoscore_cutoff = 0.95, max_kmer = 200, log_L_diff_cutoff = 0.5, psam_thres = 0.7, tmp_dir = './', name = 'RBFOX2'):
    ''' iterate '''
    
    # create logger with 'spam_application'
    logger = logging.getLogger('EM')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(
        tmp_dir,'{}_EM.log'.format(name))
    )
    fh.setLevel(logging.DEBUG)
    # add the handlers to the logger
    logger.addHandler(fh)


    # filtering dbscan result
    df = df.loc[df['ip_site'].map(len) > 0] # TODO: eliminat those clusters from DBSCAN
    df.drop(df.loc[df.iloc[:, :8].duplicated()].index, axis = 0, inplace = True) # TODO: check why there are duplicated entries
    
    # initialize psams
    
    nclus = []
    iterpsams = []
    probs = []
    seqs = []
    
    for niter in range(niters):
        logger.info('=====================iteration {}===================='.format(niter))
        if 'INCLUDE' not in df.columns:
            # if first iteration, use the best pseudoscore clusters
            
            df.loc[df['pseudoscore']>pseudoscore_cutoff, 'INCLUDE'] = True
            nclus.append(df['INCLUDE'].count())
            
        
        sub_d = df.loc[df['INCLUDE'] == True]
        logger.info('Using {} cluster to get k-mer'.format(df['INCLUDE'].count()))
        
        # calculate PSAM                
        mean, std = simulate_kmer_background(ctrl_df['seq'].tolist(), k=k, n_iter = 100, n_sample = sub_d.shape[0])
        z_score = kmer_zscore(sub_d['seq'].tolist(), mean, std, k=k)
        zcut=zscore_cutoff(z_score, plot = True, score = 'Z-score', save = os.path.join(tmp_dir,'{}_zscore_{}.png'.format(name, niter)))
        if not zcut: # no k-mer is significant
            break
            
        kmer_set = get_kmer_set(z_score, z_score_thres = zcut)
        logger.info('Total k-mer significant: {}'.format(len(kmer_set)))
        if len(kmer_set)>max_kmer:
            kmer_set = kmer_set[:max_kmer]
        
        ps=PSAMSetBuilder(kmer_set) # use top 100 7-mer to build psams
        psams = ps.make_PSAM_set()
        iterpsams.append(psams)
        # TODO set a limit to # psams
        
        # print consensus
        logger.info('PSAM found: {}'.format(' '.join([p.consensus for p in psams])))
        
        # find window
        df, use_seq = calculate_cluster_score(df, psams, psam_thres = psam_thres) # SLOW need optimize
        
        
        logger.info('PSAM window shrink to :{}'.format(' '.join(use_seq)))
        
        # learn crosslink preference
        all_p = learn_crosslink_preference(df)
        probs.append(all_p)
            
        
        
        # calculat logL
        whole_likelihood(df, all_p) # SLOW NEED OPTIMIZE
        logL_cutoff = zscore_cutoff(df['log_L'], diff_cutoff = log_L_diff_cutoff, plot = True, score = 'logL', save = os.path.join(
            tmp_dir,'{}_logL_{}.png'.format(name, niter))
        )
        #if not logL_cutoff:
        #    break # no significant crosslinking
        
        # update
        update(df, use_seq, log_L_thres= logL_cutoff, cluster_stamp = niter)
        
        
        nclus.append(df['INCLUDE'].count())
        seqs.append(use_seq)
        
    return df, nclus, iterpsams, probs, seqs

def plot_growth(nclus, save = None, title = 'cluster growth'):
    f = plt.figure()
    plt.plot(nclus, 'o-')
    plt.title(title)
    plt.xlabel('# iter')
    plt.ylabel('# clusters')

    if save:
        plt.savefig(save)
    plt.close(f)
from RBPamp.affinitylogo import plot_afflogo
def psams_together(psams, title = '', fname = ''):
    f, axes = plt.subplots(len(psams),1, figsize = (4, len(psams)*1), sharex = True)
    n = 0
    for psam, ax in zip(psams, axes):
        print('plotting psam')
        plot_afflogo(ax,psam.psam, title = '')
        if n < len(axes)-1:
            # not last subplot
             ax.set_xticklabels([])
             ax.set_xlabel('')
    plt.suptitle(title)
    f.savefig(fname)
    plt.close(f)
    
def plot_psam(iterpsams, save = None):
    for niter,psams in enumerate(iterpsams):
        print('running for {} iteration'.format(niter))
        f = psams_together(psams, title = 'iter: {}'.format(niter), fname = save+str(niter)+'.png')
        

def plot_crosslink(probs, seqs, save = None):
    f, axes = plt.subplots(len(seqs), 1, figsize = (6,len(seqs)*3))
    n=0
    print(seqs)
    for p, ax, seq in zip(probs, axes, seqs):
        labels = np.empty(p.shape, dtype = "<U1")
        zero_col = list(p.columns).index(0)
        for pid, s in enumerate(seq):
            if 'psam_{}_score'.format(pid) in p.index:
                # if not in p, then at that round no one belongs to that motif and therefore no crosslink preference
                rowid = list(p.index).index('psam_{}_score'.format(pid))
                #print(pid, zero_col, len(s), labels.shape)
                labels[rowid, zero_col:zero_col+ len(s)] = np.array(tuple(s))
        sns.heatmap(p, ax = ax, annot = labels, fmt = '')
        ax.set_title('iteration {}'.format(n))
        plt.setp(ax.get_xticklabels(), visible=False) 
        plt.setp(ax.get_yticklabels(), visible=False)
        n += 1 
    f.savefig(save+'.png')
    plt.close(f)

def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS CLIPPER FOR REGIONAL SHRINKING
        python em_refine.py --clus <pandas pickle> --ctrl <pandas pickle> --outdir <fname> --tmpdir <directory>
        
        To refine cluster boundry and poorly supported cluster using motif and crosslinking preference
        """
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012, Hsuan-lin Her 2020.
                         CLIP peakfinder by EM refinement."""
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--clus", "-i", dest="clus", help="DBSCAN pickle", type="string")

    parser.add_option("--ctrl", "-s", dest="control_clus", help="DBSCAN control pickle")
    
    parser.add_option("--outdir", "-o", dest="outdir", help="output file prefix: /home/hsher/rbfox2")
    parser.add_option("--tmpdir", "-d", dest="tmpdir", help="output file prefix: /home/hsher/rbfox2")
    
    parser.add_option("--maxkmer", dest="max_kmer", help="maximal number of kmers used to generate PSAMs", default = 200, type = "int")
    parser.add_option("--kmer", "-k", dest="k", help="length of kmer", default = 7, type = "int")
    parser.add_option("--niter", "-n", dest="niter", help="number of iteration for EM", default = 10, type = "int")
    parser.add_option("--pseudoscore", dest="pseudoscore", help = 'threshold for pseudoscore, determine initializing clusters', default = 0.95, type ="float")
    parser.add_option("--logL", dest="l", help = 'threshold for logL for crosslinking, determine cluster inclusion', default = 0.1, type ="float")
    parser.add_option("--psam", dest="psam",help = 'threshold for shrinking PSAMs, determine PSA< discrimination', default = 0.7, type ="float")

    parser.add_option("--pool", "-p", dest="pool", help="multiprocessing", default = 8, type = "int")
    
    
    
    
    return parser

if __name__ == '__main__':
    parser = option_parser()
    (options, args) = parser.parse_args()

    # load DBscan clusters and control clusters
    print('load df')
    df = pd.read_pickle(options.clus)
    name = options.clus.split('/')[-1].replace('.pickle', '')
    print('name: ', name)
    ctrl_df = pd.read_pickle(options.control_clus)

    # TODO fix dbscan to make cluster boarder and seq coherent; this is temporary
    df['start'] = df['start']-10
    df['end'] = df['end']+10

    # TODO parallelize to make this faster
    smaller_df = df.loc[(df['chrom']=='chr1') | (df['chrom']=='chr2') | (df['chrom']=='chr3')]
    print('input df shape {}'.format(smaller_df.shape))
    
    print('running em')
    result_df, nclus, iterpsams, probs, seqs = main(
        smaller_df, ctrl_df, k = options.k, niters = options.niter, 
        pseudoscore_cutoff = options.pseudoscore, max_kmer = options.max_kmer, 
        log_L_diff_cutoff = options.l, psam_thres = options.psam, tmp_dir = options.tmpdir, name = name)
    
    # save result
    print('DONE! result has {} included clusters'.format(result_df['INCLUDE'].count()))
    print('Growth: nclus:', nclus)
    result_df.to_csv(os.path.join(options.outdir, name+'.csv'))
    
    # plot cluster growth
    plot_growth(nclus, save = os.path.join(options.tmpdir, name+'_grow.png'), title = name + 'growth curve')

    # plot psams per iter
    plot_psam(iterpsams, save = os.path.join(options.tmpdir,name + '_psam.png'))

    # plot crosslinking result
    print(seqs)
    plot_crosslink(probs, seqs, save = os.path.join(options.tmpdir,name+'_crosslink'))

    print('intermediate result written to {}'.format(options.tmpdir))

