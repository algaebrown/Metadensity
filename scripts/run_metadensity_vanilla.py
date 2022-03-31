# get from metadensity from command line
import os
import pandas as pd
import metadensity as md
import warnings
import matplotlib.pyplot as plt
from optparse import OptionParser


def get_transcripts_id(eCLIP_obj, field = 'gene_id'):
    ''' return transcript IDs containing transcripts to include'''

    if transcript_ids is None:
        try:
            t = transcript.intersect(eCLIP_obj.idr, s = True, wa = True).filter().saveas()
            print('used IDR peaks containing genes')
        except Exception as e:
            print(e)
            try:
                t = transcript.intersect(eCLIP_obj.peaks['rep1'].intersect(
                    eCLIP_obj.peaks['rep2'], s = True),s = True, u = True)
                print('Fall back to: using Genes containing peaks from both replicates')
            except Exception as e:
                print(e)
                t = transcript.intersect(eCLIP_obj.peaks['rep1'],s = True, u = True)
                print('Fall back to: using Genes containing peaks from rep1 peaks')
    return [i.attrs[field] for i in t]


def main(data_row, outdir='.', single_end = True, transcript_ids = None, 
        truncation = False, bg_method = 'relative information',
        normalize = False, stat = 'mean'):
    
    
    e = eCLIP.from_series(data_row, single_end = single_end)

    if transcript_ids is None:
        transcript_ids = get_transcripts_id(e)

    if truncation:
        m = Metatruncate(e, e.name ,
                        transcript_ids = transcript_ids,
                        background_method = bg_method, 
                        normalize = normalize, 
                        sample_no = len(transcript_ids))
        m.get_density_array(use_truncation = True)
    else:
        m = Metadensity(e, e.name,
                        transcript_ids = transcript_ids,
                        background_method = bg_method, 
                        normalize = normalize,
                        sample_no = len(transcript_ids))
        m.get_density_array()
    
    print(f'Using {len(m.metagene.keys())} genes/transcripts')
    not_found = set(transcript_ids)-set(m.metagene.keys())
    
    with open(os.path.join(outdir, 'GENE_ID_NOT_FOUND.txt'),'w') as f:
        for id_ in list(not_found):
            f.write(id_+'\n')
    

    # save file  
    h5_path = os.path.join(outdir, '{}.h5'.format(e.uID))
    m.save_deepdish(h5_path)
    print(f'raw data saved to {h5_path}')
    
    # plot

    print('plotting mean density')
    
    
    for feat, featname in zip([generic_rna, protein_coding, branchpoints], ['RNA', 'CDS', 'br']):
        f = plot_mean_density([m], ymax = 0.0007, 
                            features_to_show = feat,
                            stat = stat)
        f = beautify(f)
        f.savefig(os.path.join(outdir,f'{e.uID}_{featname}.pdf'))

        f = plot_rbp_map([m], features_to_show = feat)
        f = beautify(f)
        f.savefig(os.path.join(outdir,f'{e.uID}_{featname}.map.pdf'))

    

    
def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS METADENSITY WITH ITS MOST VANILLA FUNCTION.
        For Fancy functions, please refer to the example notebooks on our website :)
        python run_metadensity.py --ip <bamfile> --input <bamfile> --feature <gff> --out <fname>
        """
    description = """Metadensity.Hsuan-lin Her 2021.
                    Pools eCLIP density across the transcriptome"""
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--csv", "-i" , dest="csv", help=".csv file containing all CLIP files", type="string")
    parser.add_option("--uid", "-u", dest="uID", help="unique ID(uid) to the CLIP in the csv you want to run on ")
    parser.add_option("--transcript_list", "-t", dest="transcript_list", 
        help="list of transcript ids to include in the calculation, if not specify, use peak-containing transcripts",
        default=None)
    parser.add_option("--out", "-o", dest="outdir", help="output path (figures and arrays)", default='.')
    parser.add_option("--single end", "-s", dest="single_end", help="Whether your CLIP is single end. Affects Metatruncate objects", action = "store_true", default = False)
    parser.add_option("--config", '-c', dest = "config", help="file path to the config file, genome coordinate specific")
    parser.add_option("--stat", dest='stat', default = 'mean', help="choose [mean,median]")
    parser.add_option("--background_method", dest='bg', help="how you want to compute IP to INPUT, choose [relative information,substraction,None]", default = 'relative information')
    parser.add_option("--normalization", dest='norm', help="whether to average the signal in a transcript", action = "store_true", default = False)
    parser.add_option("--truncation", dest='trun', help="Use truncation instead of the entire read", action = "store_true", default = False)
    return parser
    
def check_file_exist(data_path):
    ''' check if file exist for files in series'''

    isfile = data_path.apply(os.path.isfile)
    missing_fields = (isfile[~isfile]).index.tolist()
    not_path = set(missing_fields)-set(['RBP', 'uid', 'Cell line'])
    

    if len(not_path) > 0:
        warnings.warn(f'fields with file path that does not exist {not_path}' )
        return False
    return True

def parse_transcript_id_file(fname):
    ''' read file containing transcript/gene ids, store in a list'''
    with open(fname) as f:
        return [l.rstrip() for l in f.readlines()]

if __name__ == "__main__": 
    
    parser = option_parser()
    (options, args) = parser.parse_args()

    # read files for eCLIP
    datacsv = pd.read_csv(options.csv, index_col = 0)
    
    # get the exact row we want
    row = datacsv.loc[datacsv['uid']==options.uID]
    if row.shape[0]>1:
        warnings.warn(f'''
        Your Unique id (field uid) is not unique. This is VERY BAD. 
        We detected 2 rows with uid {options.uID}
        Please modify your .csv file!
        ''')
    else:
        data_path = row.iloc[0]
    print(type(data_path))
    print(data_path.apply(os.path.isfile))
    
    check_file_exist(data_path)

    # load config    
    configpath = options.config
    md.settings.from_config_file(configpath)

    from metadensity.metadensity import *
    from metadensity.plotd import *

    # read transcript ID list
    if options.transcript_list is not None:
        transcript_ids = parse_transcript_id_file(options.transcript_list)
    else:
        transcript_ids = None
    
    if options.bg == 'None':
        options.bg=None


    main(data_path, outdir=options.outdir,
     single_end = options.single_end, transcript_ids = transcript_ids, 
        truncation = options.trun, bg_method = options.bg,
        normalize = options.norm, stat = options.stat)



