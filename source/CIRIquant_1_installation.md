# Installation

**NOTES: Only python2 is supported**

## Prerequisites

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

## Install CIRIquant from source code

**Please use the latest released version from [GitHub](https://github.com/Kevinzjy/CIRIquant/releases) or [SourceForge](https://sourceforge.net/projects/ciri/files/CIRIquant/)**

Use the setup.py for CIRIquant installation (clean install using `virutalenv` is highly recommended).

```bash
# create and activate virtual env
pip install virtualenv
virtualenv venv
source ./venv/bin/activate

# Install CIRIquant and its requirement automatically
tar zxvf CIRIquant.tar.gz
cd CIRIquant
python setup.py install

# Manual installation of required pacakges is also supported
pip install -r requirements.txt
```

The package should take approximately 40 seconds to install on a normal computer.

## Install CIRIquant using pip

```
pip install CIRIquant
```
