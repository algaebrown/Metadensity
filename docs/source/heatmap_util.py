import numpy as np
import seaborn as sns
from metadensity.plotd import *
# let's just load 1 example 
base_dir = '/home/hsher/densities'

################## I KNOW THIS IS INELEGANT BUT WE STILL NEED TO SET COLORS #################
color_dict = {'5\' UTR start': 'tan',
             '5\' UTR end': 'darkgoldenrod',
              'first exon start': 'aquamarine',
              'first exon end': 'mediumaquamarine',
              'exon start': 'mediumseagreen',
              'exon end': 'seagreen',
              'last exon start': 'limegreen',
              'last exon end': 'forestgreen',
              
              'first CDS start': 'rebeccapurple',
              'first CDS end': 'mediumpurple',
              'CDS start': 'darkorchid',
              'CDS end': 'mediumorchid',
              'last CDS start': 'slateblue',
              'last CDS end': 'darkslateblue',
              'full CDS start': 'skyblue',
              'full CDS end': 'lightblue',
              
              'intron start': 'cadetblue',
              'intron end': 'darkcyan',
              
              '3\' UTR start': 'hotpink',
              '3\' UTR end': 'deeppink',
              'branchpoint start':'tomato',
              'branchpoint end':'darkorange',
              'machine predicted branchpoint start': 'lightcoral',
              'machine predicted branchpoint end': 'rosybrown',
              'polyA site start': 'grey',
              'polyA site end': 'darkgrey',
              'polyA signal start': 'wheat',
              'polyA signal end': 'goldenrod',
     }
key_to_name = {'exon': 'exon',
 'last_CDS': 'last CDS',
 'intron': 'intron',
 'full_CDS': 'full CDS',
 'first_CDS': 'first CDS',
 'last_exon': 'last exon',
 'three_prime_UTR': '3\' UTR',
 'five_prime_UTR': '5\' UTR',
 'first_exon': 'first exon',
 'CDS': 'CDS',
    'branchpoint':'branchpoint',
    'branchpoint_pred': 'machine predicted branchpoint',
    'polyAsite': 'polyA site',
              'polyAsignal': 'polyA signal'}

def read_precomputed_array(uid, base_dir = base_dir, suffix = 'densityarr'):
    ''' load deep dish h5 into data structure'''
    return dd.io.load(os.path.join(base_dir, '{}_{}.h5'.format(uid, suffix)))
def get_feature_length(denarray):
    ''' return feature length for each feature in dictionary'''
    features = list(set([d[0] for d in denarray.keys()]))
    feature_len = [denarray[f, 'left', 'rep1'].shape[1] for f in features]
    return dict(zip(features, feature_len))
def merge_two_reps(denarray, features = generic_rna, stat = np.nanmean):
    ''' merge two rep, return mean/median density'''
    all_values = []
    for f in features:
        for align in ['left', 'right']:
            
            rep1 = denarray[f, align, 'rep1']
            rep2 = denarray[f, align, 'rep2']
            
            all_den = np.concatenate([rep1,rep2], axis = 0) # concat into 1 array
            values = stat(all_den, axis = 0)
            all_values.append(values)
    return(np.concatenate(all_values))
def into_one_df(uids, features = generic_rna, base_dir = base_dir, suffix = 'densityarr', stat = np.nanmean):
    ''' given a list of uid, concat all features and data into 1 df'''
    all_vector = []
    success_uids = []
    i = 0
    for uid in uids:
        try:
            denarray = read_precomputed_array(uid, base_dir = base_dir, suffix = suffix)
            if i == 0:
                flen = get_feature_length(denarray)
            all_vector.append(merge_two_reps(denarray, features, stat = stat)) # append into results
            success_uids.append(uid)
            i+=1
        except Exception as e:
            print(e)
            print(uid)
    
    
    
    df = pd.DataFrame(np.stack(all_vector), index = success_uids)
    df.fillna(0, inplace = True)
    df.replace(np.nan, 0, inplace = True)
    
    return df , flen

def get_feat_color(features, flen):
    ''' generate color annotation for heatmap'''
    # get feature length
    colors = []
    names = []
    for f in features:
        length = flen[f]
        name = key_to_name[f]
        
        for align in ['start', 'end']:
            color_key = '{} {}'.format(name, align)
            names.append(color_key)
            colors += [color_dict[color_key]]*length
        
    return colors, names