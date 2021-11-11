# this script prepares RNAfold needed files
# needs to be run under Metadensity environment

import sys
import os
sys.path.append('/home/hsher/Metadensity/scripts')
from dataloader import *
from shapeloader import *
from metadensity.sequence import *
from metadensity.truncation import *
from pybedtools import BedTool
import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def fetch_data(chrom, start, end, strand, ip1 = None, ip2 = None, in1 = None, in2 = None, shape = None, single_end = False, read2 = True):
    ''' fetch truncation, shape and sequence for you given an interval'''

    # SHAPE
    struct = shape.fetch(chrom = chrom, start = start, end = end, strand = strand)

    # truncations
    trun1 = truncation_relative_axis(bam_fileobj=ip1, chrom = chrom, start = start, end = end, strand = strand, single_end = single_end, read2 = read2)
    trun2 = truncation_relative_axis(bam_fileobj=ip2, chrom = chrom, start = start, end = end, strand = strand, single_end = single_end, read2 = read2)
    trunin1 = truncation_relative_axis(bam_fileobj=in1, chrom = chrom, start = start, end = end, strand = strand, single_end = single_end, read2 = read2)
    all_trun = [trun1, trun2, trunin1]
    if in2:
         trunin2 = truncation_relative_axis(bam_fileobj=in2, chrom = chrom, start = start, end = end, strand = strand, single_end = single_end, read2 = read2)
         all_trun.append(trunin2)

    seq = get_interval_seq(chrom = chrom, start = start, end = end, strand = strand)

    return seq, struct, all_trun
    

def struct_to_shape_file(struct, fout = '~/scratch/shape.dat'):
    ''' format SHAPE file into shape.dat file for RNAfold '''
    df = pd.DataFrame([list(range(1,len(struct)+1)), struct]).T
    
    df.fillna(-999, inplace = True)
    df.to_csv(fout, sep = '\t', index = False, header = False)
    return df

if __name__ == '__main__':

    # read the clusters file
    cluster_file = sys.argv[1]
    if cluster_file.endswith('pickle'):
        clusters = pd.read_pickle(cluster_file)
        bed = False
    else:
        clusters = BedTool(cluster_file)
        bed = True

    tempdir = sys.argv[2]

    # get relevant information
    uid = os.path.basename(cluster_file).split('.')[0].split('_')[0]
    try:
        row = encode_data.loc[encode_data['uid']==uid]
        
    except:
        row = encode4_data.loc[encode4_data['uid']==uid]
    
    cell_line = row['Cell line'].iloc[0]
    rbp = row['RBP'].iloc[0]

    # load SHAPE data
    shape = SHAPE_data(cell_line)

    # get bam files # TODO this only work for encode3
    try:
        ip1, ip2, in1= return_fobj3(uid)
        in2 = None
        single_end = False
    except:
        ip1, ip2, in1, in2 = return_fobj4(uid)
        single_end = True

    # TODO parameter search
    window = 50
    pseudoscore_thres = 0.7
    max_data = 20

    if not bed:
        bedintervals = BedTool.from_dataframe(clusters.loc[clusters['pseudoscore']>pseudoscore_thres, ['chrom', 'start','end', 'seq', 'pseudoscore','strand']])
        
    else:
        bedintervals = clusters # bedtool
    
    for interval in bedintervals:
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        strand = interval.strand

        
        middle = math.ceil((start+end)/2)
        start = middle-window
        end = middle + window
        fname = f'{uid}_{rbp}_{cell_line}_{chrom}_{start}_{end}_{strand}'

        seq, struct, all_trun = fetch_data(chrom, start, end, strand, ip1 = ip1, ip2 = ip2, in1 = in1, in2 = in2, shape = shape, single_end = single_end)
        # check if data is too sparse
        i = 0
        if not np.all(np.isnan(struct)):

            # write file
            struct_to_shape_file(struct, fout = os.path.join(tempdir, f'{fname}.shape.dat'))
            SeqIO.write(SeqRecord(seq, id = fname), handle = os.path.join(tempdir, f'{fname}.fasta'), format = 'fasta')
            if len(all_trun) == 3:
                df = pd.DataFrame([list(seq), struct] + all_trun, index = ['seq', 'struct', 'ip1', 'ip2', 'in1'])
            else:
                df = pd.DataFrame([list(seq), struct] + all_trun, index = ['seq', 'struct', 'ip1', 'ip2', 'in1', 'in2'])

            
            df.T.to_csv(os.path.join(tempdir, f'{fname}.csv'))

            print(f'file write to {fname}')
            print('RUN conda activate rnatools')
            print(f'RNAfold --shape {fname}.shape.dat < {fname}.fasta > /home/hsher/scratch/fold.fa')
            i+=1
        
        if i > max_data:
            break




    
    
