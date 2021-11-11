import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
def score_cum_dis(df, motif, title = 'RBFOX2'):
    df['motif'] = df['seq'].str.contains(motif)

    contain_motif = df.loc[df['motif']]
    no_motif = df.loc[df['motif'] == False]
    print(no_motif.shape)
    f, ax = plt.subplots(1,1)
    
    for niter in range(10):
        smpl = no_motif.sample(contain_motif.shape[0]) # sample equal datapoints
        smpl.sort_values(by = 'pseudoscore', inplace = True)
        ax.plot(smpl['pseudoscore'].values, color = 'grey', label = 'non-motifer')

    
    contain_motif.sort_values(by = 'pseudoscore', inplace = True)    
    ax.plot(contain_motif['pseudoscore'].values, color = 'tomato', label = 'motifs')

    ax.set_title(title + '_' + motif)
    ax.set_xlabel('rank by pseudoscore')
    ax.set_ylabel('pseudoscore ')
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
#plt.plot(scores, label = 'all cluster')
#plt.plot(canon_scores, label = 'cluster with GCAUG or UCGAU')

df = pd.read_pickle('/home/hsher/encore_region_shrink/753_TIAL1_HepG2.pickle')
df_control = pd.read_pickle('/home/hsher/encore_region_shrink/753_TIAL1_HepG2.ctrl.pickle')
plt.style.use('seaborn-white')

score_cum_dis(df, 'GCAUG')

motifs = {'FUS':['GGGGG'],
        'HNRNPK': ['GCCCA'], 
        'HNRNPL': ['ACACA'],
        'KHSRP': ['UGUAU', 'UAUAU'],
        'PCBP2': ['CCCCU', 'UCCCC'],
        'EWSR1': ['GGGGG'],
          'TAF15': ['GGGGG', 'AGGGG'],
          'IGF2BP1': ['AUACA'],
          'RBFOX2': ['GCAUG', 'UGCAU'],
          'TARDBP': ['GUAUG'],
          'PUM1':['UGUAU', 'UAUAU'],
          'FUBP3':['UAUAU'],
          'PTBP3':['UUUCU'],
          'KHDRBS2':['AUAAA', 'AAUAA'],
          'SFPQ':['UGUAA'],
          'SRSF9':['AGGAG'],
         'HNRNPC': ['UUUUU', 'UUUUG'],
         'TRA2A':['GAAGA'],
         'TIA1': ['UUUUU'],
          'IGF2BP2':['ACACA'],
          'RBM22':['ACCGG'],
          'PABPN1L':['AAAAA'],
          'RBM15B':['UUUUA']}
          