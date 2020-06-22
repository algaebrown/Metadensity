# get from metadensity from command line
print('uid', 'outdir')
import sys
import os
from metadensity import eCLIP, pos_spec_bind_strength, pos_relative_entropy
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
    _ = [m.quantile_metadensity_raw(e, q = 40) for m in e.idr_metagene.values()]
    e.get_density_array(use_quantile = True)
    
    # scale 
    e.scale_density_array('max_scalar', quantile = True)
    
    print('calculate prob and entropy')
    # calculate probability
    e_prob = pos_spec_bind_strength(e, use_quantile = True, bins = 40, peak_max = 40)
    
    # get relative entropy
    #e_entropy = pos_relative_entropy(e_prob)
    
    print('saving results to {}'.format(outdir))
    
    # save aligned density
    dd.io.save(os.path.join(outdir, '{}_scaled_qdensityarr.h5'.format(e.uID)), eCLIP.scaled_density_array)
        
    # save probability
    dd.io.save(os.path.join(outdir, '{}_prob.h5'.format(e.uID)), e_prob)
    
    # save entropy
    #dd.io.save(os.path.join(outdir, '{}_entropy.h5'.format(e.uID)), e_entropy)

if __name__ == "__main__": 
    main(uid, outdir)


