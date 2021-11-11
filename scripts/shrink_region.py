from metadensity.shrink_region import *
from optparse import OptionParser
def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS CLIPPER FOR REGIONAL SHRINKING
        python shrink_region.py --ip <bamfile> --input <bamfile> --bed <bed> --out <fname>
        
        To shrink enriched regions into binding sites that we hope to contain motifs
Returns clusters (putative binding sites) and control clusters (size-matched, region-matched).
By comparing cluster to control clusters, we can calculated R values.
Clustering algorithm is DBSCAN
        """
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012, Hsuan-lin Her 2020.
                         CLIP peakfinder by shrinking enriched region."""
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--ip", "-i", dest="IPbam", help="IP bam file", type="string", metavar="FILE.bam")

    parser.add_option("--input", "-s", dest="SMInputbam", help="size match input bam file")
    
    parser.add_option("--bed", "-b", dest="enrichBed", help="bedfile output by region_call.py")
    parser.add_option("--tsv",  dest="tsv", help="tsv from Evan's pipeline", default = None)
    parser.add_option("--out", "-o", dest="outfile", help="output file prefix: /home/hsher/rbfox2")
    parser.add_option("--maxfeat", "-n", dest="nfeat", help="maximal number of feature", type = "int")
    parser.add_option("--thread", "-t", dest="thread", help="multiprocessing", default = 8, type = "int")
    parser.add_option("--eps", "-e", dest="eps", help="DBSCAN param", default = 2, type = "int")
    parser.add_option("--min_samples", "-m", dest="min_samples", help="DBSCAN param", default = 2, type = "int")
    parser.add_option("--cluster_size", "-c", dest="size", help="cluster_size", default = 30, type = "int")
    parser.add_option("--ri", "-r", dest="ri", help="use RI", action = "store_true", default = False)
    parser.add_option("--use_read1", dest="read2", help="use read1's trunction sites", action = "store_false", default = True)
    
    
    
    
    return parser

if __name__ == "__main__":
    parser = option_parser()
    (options, args) = parser.parse_args()
    # call main
    if options.tsv:
        print('Using Evan enrichment')
        clus_df, ctrl_df = main(options.IPbam, options.SMInputbam, tsv=options.tsv, n_upper = options.nfeat, 
        n_pool = options.thread, timeout = 1000, eps = options.eps, min_samples = options.min_samples, size = options.size, use_ri = options.ri, read2 = options.read2)
    else:
    
        clus_df, ctrl_df = main(options.IPbam, options.SMInputbam, enrich_bed = options.enrichBed, n_upper = options.nfeat, 
        n_pool = options.thread, timeout = 1000, eps = options.eps, min_samples = options.min_samples, size = options.size, use_ri = options.ri, read2 = options.read2)
    
    # save file
    print(clus_df.head())
    print(ctrl_df.head())

    print('saveto', options.outfile)
    clus_df.to_pickle(options.outfile + '.pickle')

    if options.ri:
        ctrl_df = ctrl_df.loc[ctrl_df['no_input']>0] ## remove those that has no crosslink
    
    print(clus_df.shape)
    print(ctrl_df.shape)
    ctrl_df.to_pickle(options.outfile + '.ctrl.pickle')
