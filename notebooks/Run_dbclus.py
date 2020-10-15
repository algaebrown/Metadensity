
def get_all_sequences(enriched_region, bam_fileobj, inputbam_fileobj, combine = False):
    ''' fetch all clusters and sequences'''
    all_clus = []
    all_control_clus = []
    all_seqs = []
    all_control_seqs = []
    for interval in enriched_region:
        clusters, control_clusters = region_cluster(interval, bam_fileobj, inputbam_fileobj, combine = combine)
        all_clus.extend(clusters)
        all_control_clus.extend(control_clusters)
        all_seqs.extend([c.seq for c in clusters])
        all_control_seqs.extend([c.seq for c in control_clusters])
    
    return all_clus, all_seqs, all_control_clus, all_control_seqs
def region_cluster(interval, bam_fileobj, inputbam_fileobj, vis = False, combine = False):
    clusters, input_sites, sites = find_cluster(interval, bam_fileobj,inputbam_fileobj, combine = combine)
    get_cluster_seq(clusters, window = 10)
    control_clusters = control_cluster(clusters, interval)
    get_cluster_seq(control_clusters, window = 10)
    if vis:
        visualize_cluster(input_sites, sites, clusters, control_clusters)
    
    return clusters, control_clusters
def get_cluster_seq(clusters, window = 10):
    ''' find cluster sequence with window'''
    seqs = []
    for c in clusters:
        
        if type(c.start)!= int or type(c.end) != int:
            print(c.start, c.end)
        s = get_interval_seq(c.chrom, c.start-window,c.end+window, c.strand)
        c.seq = s
def call_clusters(bam, inputbam, bed, outfname, n_region = 500):
    ''' given bam, input bam and the region enrich bed file'''
    # create file obj
    enriched_region = BedTool(bed)
    bam_fobj = pysam.Samfile(bam, 'rb')
    inputbam_fobj = pysam.Samfile(inputbam, 'rb')
    
    # calling clusters
    all_clus, all_seqs, all_control_clus, all_control_seqs = get_all_sequences(enriched_region[:n_region], bam_fobj, inputbam_fobj, combine = True)
    
    # save clusters
    clus_df = pd.DataFrame([[c.chrom, c.start, c.end, c.strand, c.no_ip, c.no_input, c.score, str(c.seq), c.ip_site, c.input_site] for c in all_clus],
                 columns = ['chrom', 'start', 'end', 'strand', 'no_ip', 'no_input', 'score', 'seq', 'ip_site', 'input_site'])
    clus_df.to_pickle('/home/hsher/encore_region_shrink/{}.pickle'.format(outfname))
    
    # save controls
    ctrl_clus_df = pd.DataFrame([[c.chrom, c.start, c.end, c.strand, c.no_ip, c.no_input, c.score, str(c.seq)] for c in all_control_clus],
                 columns = ['chrom', 'start', 'end', 'strand', 'no_ip', 'no_input', 'score', 'seq'])
    ctrl_clus_df.to_pickle('/home/hsher/encore_region_shrink/{}_ctrl.pickle'.format(outfname))
    
    # to bed files
    #whole_str = ''
    #for c in all_clus:
        #whole_str += c.to_bedstr()
    
    #with open('/home/hsher/encore_region_shrink/{}.bed'.format(outfname), 'w') as f:
        #f.write(whole_str)

import os
with_enrich = [f.replace('.bed', '')+'.bam' for f in os.listdir('/home/hsher/encore_region_call') if f.endswith('bed')]
subset_rbp = encode_data.loc[(encode_data['bam_0'].isin(with_enrich)) | (encode_data['bam_1'].isin(with_enrich))]
from os import path
for i in subset_rbp.index:
    inputbam = '/home/hsher/seqdata/eclip_raw/{}'.format(subset_rbp.loc[i]['bam_control'])
    uid = subset_rbp.loc[i]['uID']
    rbp = subset_rbp.loc[i]['RBP']
    if path.exists('/home/hsher/encore_region_shrink/{}_{}_{}.pickle'.format(uid, rbp, 'rep1')):
        pass
    else:
    
        # rep1
        fname = subset_rbp.loc[i]['bam_0']
        bam = '/home/hsher/seqdata/eclip_raw/{}'.format(fname)
        bed = '/home/hsher/encore_region_call/'+fname.replace('.bam', '.bed')
        try:
            call_clusters(bam, inputbam, bed, '{}_{}_{}'.format(uid, rbp, 'rep1'))
        except:
            print(bed)
    if path.exists('/home/hsher/encore_region_shrink/{}_{}_{}.pickle'.format(uid, rbp, 'rep2')):
        pass
    else:
        # rep2
        # rep1
        fname = subset_rbp.loc[i]['bam_1']
        bam = '/home/hsher/seqdata/eclip_raw/{}'.format(fname)
        bed = '/home/hsher/encore_region_call/'+fname.replace('.bam', '.bed')
    
        try:
            call_clusters(bam, inputbam, bed, '{}_{}_{}'.format(uid, rbp, 'rep2'))
        except:
            print(bed)
    print(rbp, uid, 'done')
    
    
    
