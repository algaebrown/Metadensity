# Metadensity
Visualize eCLIP/STAMP metagene density to predict RNA binding protein function

# Install
```
git clone https://github.com/algaebrown/Metadensity.git
cd Metadensity
conda create -n Metadensity
conda activate Metadensity
python setup.py install
```

# Four simple steps to create your density
1. create experimental object `eCLIP()` or `STAMP()` to hold your .bigwig/.bed/.bam files
2. select the set of transcripts you want to use using `transcript.filter(lambda x: x.attrs['gene_id'].isin(YOUR_GENE_ID_LIST)`
3. create `Metadensity()` or `Metatruncate()` object this will pool features for you. Also specify how you wannt deal with background and normalization with this object.
4. plot it using `plotd.plot_mean_density()` or `plotd.plot_rbp_map()`. Specify the set of features with `features_to_show=['five_utr', 'exon', 'three_utr']`


# Genomic Coordinate?
Default is Gencode.v33. You can generate your custom corrdinate with [this script](https://github.com/algaebrown/ClipNet/blob/master/scripts/gencode_canon_filtering.sh)

It filters for known canonical transcripts. To download it visit [UCSC table broswer](https://genome.ucsc.edu/cgi-bin/hgTables)
Download the comprehensive gene annotation gff3 [here in gencode](Comprehensive gene annotation)

Of course you can decide not to filter for known canonical transcripts.
But it will result in lots of overlapping transcript. You might see a bunch of metagene density that look the same (since they come from the same region essentially)

To point to your custom genome coordinate files visit the example notebook:[here](https://github.com/algaebrown/Metadensity/blob/master/example%20analysis/Oncogenic%20pathways.ipynb)

# Example Notebooks
[here](https://github.com/algaebrown/Metadensity/tree/master/example%20analysis)
and [here](https://github.com/algaebrown/Metadensity/blob/master/notebooks/1_Normalized%20Metadensity.ipynb)
