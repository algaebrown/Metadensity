# get from metadensity from command line

import sys
import os
from metadensity.metadensity import *
from metadensity.plotd import *
import pandas as pd
import deepdish as dd
from dataloader import encode_data, encode4_data
#import sys
#sys.path.append('/home/hsher/projects/Metadensity')
import metadensity as md

# customize feature length
md.settings.feat_len['first_CDS']=75
md.settings.feat_len['CDS']=75
md.settings.feat_len['last_CDS']=75
md.settings.feat_len['intron']=15000
# import
import pandas as pd
import matplotlib.pyplot as plt

from metadensity.metadensity import *
# make sure it's v19
transcript

def eclip_from_result_dir(result_dir, idr_path):
    ''' generate series from path'''
    # find prefix
    files = [f for f in os.listdir(result_dir) if f.endswith('.rep1_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.blacklist-removed.narrowPeak')]

    prefixes = [f.split('.umi')[0] for f in files] # TRNAU1AP_2.0.rep1_CLIP, TRNAU1AP_2.0.rep2_CLIP
    reps = [0 if 'rep1' in p else 1 for p in prefixes]

    data = pd.Series()

    # uid and name
    data['uid']=prefixes[0].split('.')[0]
    data['RBP']=prefixes[0].split('.')[0]
    
    for prefix, rep in zip(prefixes, reps):
        
        #bams
        # IP
        data[f'bam_{rep}'] = os.path.join(result_dir, f'{prefix}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam')

        # INPUT
        input_prefix=prefix.replace('CLIP', 'INPUT')
        data[f'bam_control_{rep}'] = os.path.join(result_dir, f'{input_prefix}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam')

        for strand, name in zip(['neg', 'pos'], ['minus', 'plus']):
            data[f'{name}_{rep}'] = os.path.join(result_dir, f'{prefix}.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.{strand}.bw')
            data[f'{name}_control_{rep}'] = os.path.join(result_dir, f'{input_prefix}.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.{strand}.bw')

        # individual peaks
        data[f'bed_{rep}']=os.path.join(result_dir, f'{prefix}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.blacklist-removed.bed')

        # IDR path
        data['idr']=idr_path
    
    # check error
    print(data)

    print(data.apply(os.path.isfile))

    return data

def main(data_row, single_end = True):
    
    
    e = eCLIP.from_series(data_row, single_end = single_end)
    try:
        t = transcript.intersect(e.idr, s = True, wa = True).filter(lambda x: x.attrs['transcript_type']=='protein_coding').saveas()
    except:
        t = transcript.intersect(e.peaks['rep1'].intersect(e.peaks['rep2'], s = True),s = True, u = True)

    # get aligned density
    mden = Metadensity(e,  e.name+'_idr',transcripts = t, background_method = 'subtract', normalize = True)
    mden.get_density_array()
    
    # get aligned density, relative information
    mden_ri = Metadensity(e, e.name+'_idr',background_method = 'relative information', normalize = False, metagenes = mden.metagene.copy())
    mden_ri.get_density_array()
    
    # get truncation
    mtru = Metatruncate(e, e.name+'_idr',background_method = 'subtract', normalize = True, metagenes = mden.metagene.copy())
    mtru.get_density_array(use_truncation = True)
    
       
    
    # plotting
    
    rna = ['first_exon', 'exon', 'intron', 'last_exon']
    protein_coding = ['five_prime_UTR', 'first_CDS', 'CDS', 'last_CDS', 'three_prime_UTR']
    branchpoints = ['branchpoint_pred', 'branchpoint']
    
    
    print('plotting mean density')
    
    for obj, objname in zip([mtru, mden, mden_ri], ['trun', 'den', 'denri']):
        for feat, featname in zip([generic_rna, protein_coding, branchpoints], ['RNA', 'CDS', 'br']):


            f = plot_mean_density([obj], ymax = 0.0007, features_to_show = feat)
            f.savefig(os.path.join(outdir,f'{e.uID}_{featname}_{objname}.pdf'), dpi = 300)

            f = plot_rbp_map([obj], features_to_show = feat)
            f.savefig(os.path.join(outdir,f'{e.uID}_{featname}_{objname}.map.pdf'), dpi = 300)

    print('saving results to {}'.format(outdir))
    
    # save aligned density
    mden.save_deepdish(os.path.join(outdir, '{}_densityarr.h5'.format(e.uID)))
    mden_ri.save_deepdish(os.path.join(outdir, '{}_ridensityarr.h5'.format(e.uID)))
    mtru.save_deepdish(os.path.join(outdir, '{}_truncatearr.h5'.format(e.uID)))

    

    
    

if __name__ == "__main__": 
    
    result_path = sys.argv[1]
    idr_path = sys.argv[2]
    outdir = sys.argv[3]
    data_row = eclip_from_result_dir(result_path, idr_path)
    main(data_row)
    
    
     # out directory
    # datacsv = pd.read_csv(sys.argv[1], index_col = 0)
    # outdir = sys.argv[2]

    # for index,data_row in datacsv.iterrows():
        
    #     try:
    #         main(data_row)
    #     except:
    #         pass
    # get uID, points to all the files and ids



