import sys
#sys.path.append('/home/hsher/projects/Metadensity')
from metadensity.metadensity import *
from metadensity.plotd import *
from metadensity.pos_enrich import *
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append('/home/hsher/Metadensity/scripts')
from dataloader import *
import seaborn as sns
from scipy.ndimage import gaussian_filter
import numpy as np
def relative_entropy(metagene, rbp, pseudocount = 0.001, sigma = 5):
    csln = metagene.sites[rbp]
    
    ip1 = csln['rep1']
    ip2 = csln['rep2']
    ipdist = (ip1+ip2)+pseudocount
    ipdist = gaussian_filter(ipdist/ipdist.sum(), sigma = sigma)
    
    if 'ctrl' in csln.keys():
        in_ = csln['ctrl']
    else:
        in1 = csln['ctrl1']
        in2 = csln['ctrl2']
        in_ = in1+in2
    
    indist = in_+pseudocount
    indist = gaussian_filter(indist/indist.sum(), sigma = sigma)
    trun = ipdist*np.log(ipdist/indist) ### relative entropy
    
    return trun

def get_metagene_sites(metagene, sigma = 5):
    all_data = []
    uids = []
    for rbp in metagene.sites.keys():
        try:
            #data = gaussian_filter((metagene.value[rbp]['rep1'] + metagene.value[rbp]['rep2'])/2, sigma = sigma)
            data = relative_entropy(metagene, rbp, sigma = sigma)
            
            data[np.where(data < 0)]=0
            all_data.append(data)
            uids.append(rbp)
        except:
            print(rbp)
    f,ax = plt.subplots(figsize = (10,4))
    all_data = pd.DataFrame(np.stack(all_data), index = uids)
    #sns.heatmap(all_data, vmin = 0,vmax = all_data.max().max()/10, cmap = 'pink_r', ax = ax)
    return all_data
    


if __name__=='__main__':
    geneid = sys.argv[1]
    outdir = sys.argv[2]

    # build metagene for particular transcript
    ids = [geneid]
    ms = Build_many_metagene(ids, transcript_type = None) 

    # build all eCLIP object
    all_eclips = []
    for index, row in encode_data.iterrows():
        try:
            if row['uid'] in encode_data['uid'].tolist():
                all_eclips.append(eCLIP.from_series(row, single_end = False))
            else:
                all_eclips.append(eCLIP.from_series(row, single_end = True))
        except Exception as e:
            print(e)
    for index, row in encode4_data.iterrows():
        try:
            if row['uid'] in encode_data['uid'].tolist():
                all_eclips.append(eCLIP.from_series(row, single_end = False))
            else:
                all_eclips.append(eCLIP.from_series(row, single_end = True))
        except Exception as e:
            print(e)
    
    # get all values
    for i in ids:
        _ = [ms[i].get_value(e, background_method = 'relative information', normalize = False, truncate = True) for e in all_eclips]

    # get all get_metagene_sites
    data = get_metagene_sites(ms[geneid])

    data.to_csv(os.path.join(outdir, geneid))

    

    
