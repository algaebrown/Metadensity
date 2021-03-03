# get from metadensity from command line
import sys
import os
from metadensity.metadensity import *
from metadensity.plotd import *
import pandas as pd
import deepdish as dd
from dataloader import encode_data, encode4_data

def main(uID, outdir):
    
        
    if uID[0] == '4' and len(uID)==4:
        # encode 4
        data_row = encode4_data.loc[encode4_data['uid']==uID].iloc[0]
        single_end = True
    else:
        # encode 3
        data_row = encode_data.loc[encode_data['uid'] == uID].iloc[0]
        single_end = False
    
    print('=============> getting eCLIP object from these files <=====================:')
    print(data_row)
    # create eCLIP object and get metadensity
    
    e = eCLIP.from_series(data_row, single_end = single_end)
    
    t = transcript.intersect(e.idr, s = True, wa = True).filter(lambda x: x.attrs['transcript_type']=='protein_coding').saveas()

    # get aligned density
    mden = Metadensity(e,  e.name+'_idr',transcripts = t, background_method = 'subtract', normalize = True)
    mden.get_density_array()
    
    # get aligned density, relative information
    mden_ri = Metadensity(e, e.name+'_idr',background_method = 'relative information', normalize = False, metagenes = mden.metagene.copy())
    mden_ri.get_density_array()
    
    # get truncation
    mtru = Metatruncate(e, e.name+'_idr',background_method = 'subtract', normalize = True, metagenes = mden.metagene.copy())
    mtru.get_density_array(use_truncation = True)
    
       
    print('saving results to {}'.format(outdir))
    
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_densityarr.h5'.format(e.uID)), mden.density_array)
    
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_ridensityarr.h5'.format(e.uID)), mden_ri.density_array)
        
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_truncatearr.h5'.format(e.uID)), mtru.truncate_array)

    # plotting
    
    rna = ['first_exon', 'exon', 'intron', 'last_exon']
    protein_coding = ['five_prime_UTR', 'first_CDS', 'CDS', 'last_CDS', 'three_prime_UTR']
    
    
    print('plotting mean density')
    
    f = plot_mean_density([mtru], ymax = 0.0007, features_to_show = rna)
    f.savefig(os.path.join(outdir,'{}_{}_{}_trun.pdf'.format(e.uID, 'RNA', 'mean')), dpi = 300)

    f = plot_mean_density([mden], ymax = 0.0007, features_to_show = rna)
    f.savefig(os.path.join(outdir,'{}_{}_{}_den.pdf'.format(e.uID, 'RNA', 'mean')), dpi = 300)
    
    f = plot_mean_density([mtru], ymax = 0.0007, features_to_show = protein_coding)
    f.savefig(os.path.join(outdir, '{}_{}_{}_trun.pdf'.format(e.uID, 'CDS', 'mean')), dpi = 300)

    f = plot_mean_density([mtru], ymax = 0.0007, features_to_show = protein_coding)
    f.savefig(os.path.join(outdir,'{}_{}_{}_den.pdf'.format(e.uID, 'CDS', 'mean')), dpi = 300)
    
    print('plotting map density')
    
    f = plot_rbp_map([mtru], features_to_show = protein_coding)
    f.savefig(os.path.join(outdir,'{}_{}_{}_trun.pdf'.format(e.uID, 'CDS', 'map')), dpi = 300)
    
    f = plot_rbp_map([mtru], features_to_show = rna)
    f.savefig(os.path.join(outdir,'{}_{}_{}_trun.pdf'.format(e.uID, 'RNA', 'map')), dpi = 300)

if __name__ == "__main__": 
    

    uid = sys.argv[1]
    # out directory
    outdir = sys.argv[2]
    main(uid, outdir)
    # get uID, points to all the files and ids



