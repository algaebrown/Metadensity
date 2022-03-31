from setuptools import setup
from setuptools import find_packages


long_description = "Metadensity - Metagene plot using CLIP-seq and STAMP, for RNAs"

setup(
    name = "metadensity",
    long_description = long_description,
    version = "0.0.1",
    packages = find_packages(),
    
    package_dir = {'metadensity': 'metadensity'},
    package_data = {
        'metadensity' : ['data/*/*']
        },

    install_requires = ['setuptools', 
                        'pysam >= 0.15.4',
                        'numpy >= 1.18.1 ',
                        'scipy >= 1.4.1',
                        'matplotlib >= 3.1.3',
                        'pybedtools >= 0.8.1',
                        'pandas >= 1.0.2',
                        'pybigwig,
                        'seaborn,
                        'biopython',
                        'scikit-learn',
                        'deepdish'
                        ],
      
    setup_requires = ["setuptools_git >= 0.3",],
    
    #scripts=['scripts/region_call.py', 'scripts/shrink_region.py'],
     

    #metadata for upload to PyPI
    author = "Hsuan-lin Her",
    author_email = "hsher@ucsd.edu",
    description = "A set of scripts visualizing eCLIP data",
    license = "GPL2",
    keywords = "CLIP-seq, metagene, bioinformatics",
    url = "https://github.com/algaebrown/Metadensity",
    
    #Other stuff I feel like including here
    #include_package_data = True
    #zip_safe = True #True I think
)