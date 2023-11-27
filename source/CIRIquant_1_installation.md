# Installation

**NOTES: Only python2 is supported**

## Prerequisites

```
Softwares:
- bwa
- hisat2
- stringtie
- samtools >= 1.9 (*older version of samtools may use deprecated parameters in `sort` and `index` commands*)

Python packages:
- PyYAML
- argparse
- pysam
- numpy
- scipy
- scikit-learn
```

## Use the released version (Recommended)

**Download the latest released version of CIRIquant from [GitHub](https://github.com/Kevinzjy/CIRIquant/releases)**

The released package is a packed conda environment including all dependencie, make sure you have installed anaconda in your environment

```bash
# Download packed package
wget https://github.com/bioinfo-biols/CIRIquant/releases/download/v1.1.3/CIRIquant_v1.1.3.tar.gz
mkdir -p CIRIquant_env
tar zxvf CIRIquant_v1.1.3.tar.gz -C CIRIquant_env

# Configuration environments (required)
conda activate ./CIRIquant_env
cd CIRIquant_env
make
```

Activate CIRIquant environment and test

```bash
conda activate ./CIRIquant_env
which CIRIquant
```

## Install CIRIquant and dependencies using conda (Not recommended)

Save the following content to a file called environment.yml:

```
name: CIRI
channels:
- defaults
- bioconda
- conda-forge
dependencies:
- bioconda::bwa=0.7.17
- bioconda::hisat2=2.2.0
- bioconda::stringtie=2.1.1
- bioconda::samtools>=1.10
- bioconda::bioconductor-edger=3.28.0
- bioconda::bioconductor-limma=3.42.0
- conda-forge::r-statmod=1.4.35
- conda-forge::r-base=3.6
- conda-forge::r-optparse=1.6.6
- python=2.7.15
- pip=20.0.2
- perl=5.26.2
- pip:
  - CIRIquant>=1.1.2
  - numexpr==2.6.9
  - numpy==1.16.4
  - pysam==0.15.2
  - PyYAML==5.4
  - scikit-learn==0.20.3
  - scipy==1.2.2
  - argparse>=1.2.1
```

After you have saved the file just run:
```
# this installs the dependencies specified in the yml file
conda env create -f environment.yml
# this activates the conda environment
conda activate CIRI

# this will return the path bwa, hisat2, stringtie or samtools are installed to
# these paths need to be specified in the CIRI configuration file when running the tool
which bwa
which hisat2
which stringtie
which samtools
```

## Install CIRIquant using pip

```
pip install CIRIquant
```

## Install CIRIquant from source code

Use the setup.py for CIRIquant installation (clean install using `virutalenv` is highly recommended).

```bash
# create and activate virtual env
pip install virtualenv
virtualenv -p /path/to/your/python2/executable venv
source ./venv/bin/activate

# Install CIRIquant and its requirement automatically
tar zxvf CIRIquant.tar.gz
cd CIRIquant
python setup.py install

# Manual installation of required pacakges is also supported
pip install -r requirements.txt
```

The package should take approximately 40 seconds to install on a normal computer.
