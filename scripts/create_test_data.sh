module load samtools
samtools view -h /projects/ps-yeolab3/encode/analysis/encode_GRCh38_v1/228_01_SF3B4.merged.r2.bam "chr22:10736171-50783662" | samtools view -Sb -> /home/hsher/Metadensity/test_data/processed_bam/SF3B4_CLIP.r2.bam
samtools view -h /projects/ps-yeolab3/encode/analysis/encode_GRCh38_v1/228_INPUT_TCCGGAGA-AGGCGAAG_L001_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam "chr22:10736171-50783662" | samtools view -Sb -> /home/hsher/Metadensity/test_data/processed_bam/SF3B4_INPUT.r2.bam

