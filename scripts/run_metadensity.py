# get from metadensity from command line
print('uid', 'outdir')
import sys
import os
sys.path.append('/home/hsher/projects/Metadensity')
from metadensity.metadensity import *
import pandas as pd
import deepdish as dd

# get uID, points to all the files and ids

uid = sys.argv[1]
# out directory
outdir = sys.argv[2]


# load IDs
encode_data = pd.read_pickle('~/projects/eclip_encode_id.pickle')
eclip_bam = pd.read_csv('/home/hsher/projects/RBP_annot/ENCODE_FINAL_ANNOTATIONS.uidsonly.txt.manifesthg38.txt', sep = '\t', header= 0)

# join data
encode_data = pd.merge(eclip_bam[['uID', 'RBP', 'Cell line']], encode_data, left_on = ['RBP', 'Cell line'], right_on = ['RBP', 'cell_line'])

def main(uID, outdir):
    data_row = encode_data.loc[encode_data['uID'] == uID]
    
    print('getting eCLIP object')
    # create eCLIP object and get metadensity
    e = eCLIP()
    e.RBP_centric_approach(data_row)
    
    # get aligned density
    mden = Metadensity(e, e.idr_transcript, e.name+'_idr',background_method = 'subtract', normalize = True)
    mden.get_density_array()
    
    # get aligned density, relative information
    mden_ri = Metadensity(e, e.idr_transcript, e.name+'_idr',background_method = 'relative information', normalize = False, metagenes = mden.metagene.copy())
    mden_ri.get_density_array()
    
    # get truncation
    mtru = Metadensity(e, e.idr_transcript, e.name+'_idr',background_method = 'subtract', normalize = True, metagenes = mden.metagene.copy())
    mtru.get_density_array()
    
       
    print('saving results to {}'.format(outdir))
    
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_densityarr.h5'.format(e.uID)), mden.density_array)
    
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_ridensityarr.h5'.format(e.uID)), mden_ri.density_array)
        
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_truncatearr.h5'.format(e.uID)), mtru.truncate_array)

if __name__ == "__main__": 
    main(uid, outdir)


