# Metadensity
Visualize eCLIP/STAMP metagene density to predict RNA binding protein function
[Read the Doc] https://metadensity.readthedocs.io/en/latest/

# Install
The environment is available at [DockerHub](https://hub.docker.com/repository/docker/algaebrown/metadensity)
```
git clone https://github.com/algaebrown/Metadensity.git
cd Metadensity
# build your own environment!
conda create -n Metadensity --file environment.yaml
conda activate Metadensity
# copy genome coordinate
cd Metadensity
pip install -e .
```
# Step 1: download data and setup paths to annotations in `config/`
Metadensity requires several annotations to work. You need to point those files in `config/*.ini`. see `config/hg38.ini` as an example.
These information are genome coordinate, species dependent. So you can keep each species with a seperate `.ini`.

|                  | Description                                                                              | Link to download                                                                            | Essential to run |
|------------------|------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|------------------|
| GENOME_FA        | fasta for the entire genome sequence                                                     | https://www.ncbi.nlm.nih.gov/genome/guide/human/                                            | YES              |
| GENCODE          | gff3 annotation of exon, intron, gene, transcripts etc                                   | https://www.gencodegenes.org/human/                                                         | YES              |
| BRANCHPOINT      | branchpoint annotation                                                                   | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4315302/                                       | NO               |
| BRANCHPOINT_PRED | machine learning predicted branchpoints                                                  | https://pubmed.ncbi.nlm.nih.gov/29092009/                                                   | NO               |
| POLYA            | annotation of polyA sites and signals from polyASite                                     | https://www.google.com/search?client=firefox-b-1-d&q=polyA.atlas.clusters.2.0.GRCh38.96.bed | NO               |
| MIRNA            | annotation of microRNA                                                                   | https://www.mirbase.org/                                                                    | NO               |
| SNORNA           | annotation of snoRNA                                                                     | http://scottgroup.med.usherbrooke.ca/snoDB/                                                 | NO               |
| LNCRNA           | annotation of lncRNA                                                                     | https://lncipedia.org/                                                                      | NO               |
| TRANSCRIPT       | gff3 annotation from GENCODE, containing only "transcript"                               |                                                                                             | YES              |
| FEATURE          | gff3 annotation from GENCODE, containing only "exon", "CDS", "UTR" and created "introns" |                                                                                             | YES              |
| DATADIR          | parsed information of above, python .pickle file                                         | https://www.dropbox.com/sh/hoya37n9pmuqd4l/AABBSpcpjFYIUFWMdIRuJtU4a?dl=0                   |                  |



# Step 2: Command line usage: the most vanilla functions
`python scripts/run_metadensity_vanilla.py -h`

```
Options:
  -h, --help            show this help message and exit
  -i CSV, --csv=CSV     .csv file containing all CLIP files
  -u UID, --uid=UID     unique ID(uid) to the CLIP in the csv you want to run
                        on
  -t TRANSCRIPT_LIST, --transcript_list=TRANSCRIPT_LIST
                        list of transcript ids to include in the calculation,
                        if not specify, use peak-containing transcripts
  -o OUTDIR, --out=OUTDIR
                        output path (figures and arrays)
  -s, --single end      Whether your CLIP is single end. Affects Metatruncate
                        objects
  -c CONFIG, --config=CONFIG
                        file path to the config file, genome coordinate
                        specific
  --stat=STAT           choose [mean,median]
  --background_method=BG
                        how you want to compute IP to INPUT, choose [relative
                        information,substraction,None]
  --normalization       whether to average the signal in a transcript
  --truncation          Use truncation instead of the entire read
```

This will run the most vanilla functions of Metadensity

Here we provide some [test data](https://www.dropbox.com/s/cgkeuqr0cjif558/test_data.tar.gz?dl=0) to run the script.

```
cd Metadensity
# download the data
wget https://www.dropbox.com/s/cgkeuqr0cjif558/test_data.tar.gz

# uncompress
tar -xvzf test_data.tar.gz

# modify paths is menifest.csv to correspond to your directory
nano test_data/menifest.csv

# run
python scripts/run_metadensity_vanilla.py -i test_data/menifest.csv -u SF3B4_test -o test_data --config=config/hg38.ini
```


# Step 3: Advanced usage in jupyter notebook
[checkout example jupyter notebooks](https://github.com/algaebrown/Metadensity/tree/master/docs/source)