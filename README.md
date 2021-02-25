# Metadensity
Visualize eCLIP/STAMP metagene density to predict RNA binding protein function

# Install
```
git clone https://github.com/algaebrown/Metadensity.git
cd Metadensity
# build your own environment!
conda create -n Metadensity
conda activate Metadensity
# copy genome coordinate
cd metadensity
cp /home/hsher/Metadensity/metadensity/data . # only work in TSCC
python setup.py install # install
```

# Four simple steps to create your density
1. create experimental object `eCLIP()` or `STAMP()` to hold your .bigwig/.bed/.bam files
2. select the set of transcripts you want to use using `transcript.filter(lambda x: x.attrs['gene_id'].isin(YOUR_GENE_ID_LIST)`
3. create `Metadensity()` or `Metatruncate()` object this will pool features for you. Also specify how you wannt deal with background and normalization with this object.
4. plot it using `plotd.plot_mean_density()` or `plotd.plot_rbp_map()`. Specify the set of features with `features_to_show=['five_utr', 'exon', 'three_utr']`


# Genomic Coordinate?
Currently the genome coordinates are available at `ls /home/hsher/Metadensity/metadensity/data` with only hg19, hg38 and mm10.

# Example Notebooks
[here](https://github.com/algaebrown/Metadensity/tree/master/example%20analysis)

