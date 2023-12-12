Installation
=======================================

Install from source
--------------------------------------
.. code-block:: bash
   git clone https://github.com/algaebrown/Metadensity.git
   cd Metadensity
   # build your own environment!
   conda env create -n Metadensity --file environment.yaml
   conda activate Metadensity
   # copy genome coordinate
   cd Metadensity
   pip install -e .

Build annoatations
--------------------------------------
1. Step 1: Download annoatations
.. code-block:: bash
   # on bash, download gencode
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.primary_assembly.annotation.gff3.gz
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz

   gunzip gencode.v40.primary_assembly.annotation.gff3.gz

2. compile config.ini file (Using gencode 40 as an example)
Example file is in `here <https://github.com/algaebrown/Metadensity/blob/master/config/hg38.ini>`

This is a file that stores all the annoatations
.. code-block:: bash
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.primary_assembly.annotation.gff3.gz
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz

   gunzip gencode.v40.primary_assembly.annotation.gff3.gz

Fill in GENOME_FA, GENCODE as the files that you just downloaded.
the fields in OTHER are optional. See here for download `link <https://github.com/algaebrown/Metadensity/tree/master#step-1-download-data-and-setup-paths-to-annotations-in-config>`

TRANSCRIPT and FEATURE is generated from the following steps. Simply replace .gff with .transcript.gff and combined.sorted.gff3
Lastly, DATADIR is also generated from the steps below. Simply put any directory that can hold lots of data.

3. Generate DATADIR
use `this notebook <https://github.com/algaebrown/Metadensity/blob/master/docs/source/parse_gencode_coords_into_data.ipynb>` to generate.
This notebooks takes CONFIG.ini and parse all the annotations into a dictionary.

4. Generate TRANSCRIPT and FEATURE from GENCODE and UCSC canonical transcripts
Download canonical transcripts `here <http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=725869661_wEKn8bV7KmsR6WJ9W6jIT45len1r&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=knownGene&hgta_table=knownCanonical&hgta_regionType=genome&position=chr1%3A11%2C102%2C837-11%2C267%2C747&hgta_outputType=primaryTable&hgta_outFileName=>`.
Use `this script <https://github.com/algaebrown/Metadensity/blob/master/scripts/gencode_canon_filtering.sh>` to generate TRANSCRIPT and FEATURE files.
.. code-block::
   module load bedtools
   bash gencode_canon_filtering.sh gencode.v40.primary_assembly.annotation.gff3 knownCanonical.txt

Run the test notebook
--------------------------------------
Run this `test notebook <https://github.com/algaebrown/Metadensity/blob/master/docs/source/1_Example_on_test_data.ipynb>`.