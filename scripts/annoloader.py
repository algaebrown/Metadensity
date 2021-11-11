# Load RBP functional annoatation
import pandas as pd
import numpy as np
from dataloader import *

rbp_class = pd.read_excel('~/projects/RBP_annot/rbp_class.xlsx', skiprows = 1)
rbp_class.set_index('Unnamed: 0', inplace = True)
rbp_class.drop('Unnamed: 1', axis = 1, inplace = True)

rbp_complex = pd.read_csv('~/projects/RBP_annot/RBP_complex_annotation.csv', index_col = 0)
go_cc = pd.read_csv('~/projects/RBP_annot/go_cc.csv', index_col = 0)
go_bp = pd.read_csv('~/projects/RBP_annot/go_cc.csv', index_col = 0)
go_mf = pd.read_csv('~/projects/RBP_annot/go_mf.csv', index_col = 0)
interpro_domains = pd.read_csv('~/projects/RBP_annot/inter_pro_domain.csv', index_col = 0)
interopro_repeats = pd.read_csv('~/projects/RBP_annot/inter_pro_repeats.csv', index_col = 0)

def annotation_to_color(df, true_color = 'dimgrey', false_color = 'linen', missing_color = 'white'):
       df = df.replace(1, True).replace(0, False)
       return df.replace(True, true_color).replace(False, false_color).replace(np.nan, missing_color)
# RBP class with color

subset = ['Splicing regulation', 'Spliceosome','3\' end processing', 
       'Ribosome & basic translation', 'RNA stability & decay',
       'microRNA processing', 'Translation regulation']
rbp_class_color = rbp_class.replace(1, 'dimgrey').replace(0, 'linen').replace(np.nan, 'linen')
rbp_complex_color = rbp_complex.replace(True, 'dimgrey').replace(False, 'linen').replace(np.nan, 'linen').drop(np.nan)

# for those that are yet un-annotated, make them white
un_anno = master_df['RBP'][~master_df['RBP'].isin(rbp_class.index)].unique()
rbp_class_color=rbp_class_color.append(pd.DataFrame(index = un_anno))
rbp_class_color.loc[un_anno, :] = 'white'