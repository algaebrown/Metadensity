from metadensity.truncation import *
from metadensity.sequence import *
import matplotlib.pyplot as plt
from collections import Counter
import sys
from dataloader import *
from shapeloader import *
from multiprocessing import Pool

def get_datapointer_from_full_prefix(prefix, uid, RBP, Cell_line):
    ''' given prefix in a dataframe, return another dataframe with full path'''
    datapointer = pd.Series()
    datapointer['uid'] = uid
    datapointer['RBP'] = RBP
    datapointer['Cell line'] = Cell_line
    
    
    datapointer['bam_0'] = prefix+'_CLIP1_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam'
    datapointer['bam_1'] = prefix+'_CLIP2_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam'
    datapointer['bam_control_0'] = prefix+'_CLIP1_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam'
    datapointer['bam_control_1'] = prefix+'_CLIP2_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam'
    
    # biuild bigwig filename
    datapointer['plus_0'] = prefix + '_CLIP1_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.pos.bw'
    datapointer['plus_1'] = prefix + '_CLIP2_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.pos.bw'
    datapointer['plus_control_0'] = prefix + '_CLIP1_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.pos.bw'
    datapointer['plus_control_1'] = prefix + '_CLIP2_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.pos.bw'
    datapointer['minus_0'] = prefix + '_CLIP1_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.neg.bw'
    datapointer['minus_1'] = prefix + '_CLIP2_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.neg.bw'
    datapointer['minus_control_0'] = prefix + '_CLIP1_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.neg.bw'
    datapointer['minus_control_1'] = prefix + '_CLIP2_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.neg.bw'
    
    # build individual peak file name
    datapointer['bed_0'] = prefix + '_CLIP1_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.blacklist-removed.bed'
    datapointer['bed_1'] = prefix + '_CLIP2_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.blacklist-removed.bed'
    
    # IDR
    # build IDR peak path
    #idr_base = '/projects/ps-yeolab5/encore/processing/encore_master_IDR_hg38/results_20210303/'
    #datapointer['idr']=idr_base + 'encore_master_'+ datapointer['Batch'] + '_IDR/results/' + datapointer['uid'] + '_CLIP1_rep1.vs.' + datapointer['uid'] + '_CLIP2_rep2.bed'
    
    return datapointer
# fetch eclip mismatch, trunction and indel. These are indication of a crosslinking (truncation is the most likely event)
def align_eclip_and_shape(p, datapointer):
    ''' fetch eCLIP crosslinking events for IP/INPUT for BedTool interval p
    Also fetch icSHAPE reacticity for the same region p
    join them together as a dataframe
    '''
    bam1 = pysam.AlignmentFile(datapointer['bam_0'])
    bam2 = pysam.AlignmentFile(datapointer['bam_1'])
    bamin1 = pysam.AlignmentFile(datapointer['bam_control_0'])
    bamin2 = pysam.AlignmentFile(datapointer['bam_control_1'])

    all_prof = []
    # get the eclip stuffs, the read is reverse complement to the reference
    for eclip_bam, name in zip([bam1, bam2, bamin1, bamin2], ['IP1', 'IP2', 'IN1', 'IN2']):
        eclip_profile = strand_specific_pileup(eclip_bam, chrom = p.chrom, start = p.start, end = p.end, strand = p.strand)
        eclip_profile.columns = [c+'_'+name for c in eclip_profile.columns]
        all_prof.append(eclip_profile)
    
    # join the 2 reps and the inputs together
    eclip_profile = pd.concat(all_prof, axis = 1)

    if p.strand == '+':
        reactivity=icshape_plus.values(p.chrom, p.start, p.end)
    else:
        reactivity=icshape_minus.values(p.chrom, p.start, p.end)
    
    # organize into series
    react_col = pd.Series(reactivity, index=np.arange(p.start, p.end))

    # put into the same dataframe
    eclip_profile['shape'] = react_col
    
    # fetch reference genome sequence
    seq = getRNAsequence(p)
    pos = np.arange(p.start, p.end)
    if p.strand == '-':
        pos = pos[::-1]
    
    eclip_profile['ref_seq'] = pd.Series(seq, index = pos)
    eclip_profile.to_csv('/home/hsher/seqdata/adar_clip_shape/{}.csv'.format(p.attrs['ID']))
    print(p.attrs)
    return True
    

if __name__=='__main__':
    # adar_prefix = '/projects/ps-yeolab5/encore/processing/batch23_20210302/results/encode4_batch23.4137'
    # load CLIP data
    adar_prefix = '/projects/ps-yeolab5/encore/processing/batch23_20210302/results/encode4_batch23.4137'
    datapointer = get_datapointer_from_full_prefix(adar_prefix, '4137', 'ADAR', 'K562')
    

    import pyBigWig
    # load icSHAPE data
    icshape_plus = pyBigWig.open('/home/hsher/icshape_data_zhang/icSHAPE/K562-plus.bw')
    icshape_minus = pyBigWig.open('/home/hsher/icshape_data_zhang/icSHAPE/K562-minus.bw')

    # load transcripts
    transcripts = BedTool('/home/hsher/gencode_coords/gencode.v33.transcript.gff3')


    # setup multiprocessing
    n_pool = 8
    timeout = 1000
    pool = Pool(int(n_pool)) #interval, ip_bamfile, input_bamfile, ip_total, smi_total, pval_thres = 0.05
    tasks = [
            (
            interval,
            datapointer
            
            )
             for interval in transcripts
             ] # pseudocount
    print('calling for {} regions:'.format(len(tasks)))
    jobs = [pool.apply_async(align_eclip_and_shape, task) for task in tasks]
    for job in jobs:
        try:
            
            
            # if return anything
            df = job.get(timeout=timeout)
        except:
            pass


    
    pool.close()



