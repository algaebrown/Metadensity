import sys
import os

gse=sys.argv[1]
target_dir=sys.argv[2]

gsennn=gse[:-3]+'nnn'

os.system(f'wget \'ftp://ftp.ncbi.nlm.nih.gov/geo/series/{gsennn}/{gse}/suppl/{gse}_RAW.tar\' -P {target_dir}')

