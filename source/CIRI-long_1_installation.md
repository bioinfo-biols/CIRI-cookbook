# Installation

## Dependency

- `gcc 4.8+` or `clang 3.4+` and `cmake 3.2+` is needed
- **Only python3 is supported**
- CIRI-long requires pysam lib, which need executable binary and header of zlib, bzip2, xz, please refer to documentation of pysam for installation instructions
- all python dependencies are listed in requirements.txt
- `samtools` version `1.9` or higher
- python binding of minimap2 

## Install CIRI-long from source code

Installation under `virtualenv` is highly recommended. You can simply clone the whole repository, then use `make` to start a complete installation

```bash
git clone https://github.com/bioinfo-biols/CIRI-long.git CIRI-long
cd CIRI-long

# Create virtual environment
python3 -m venv venv

# Activate virtualenv
source ./venv/bin/activate

# Install CIRI-long
make

# Test for installation
make test
```

The package should take less than 1 min to install on a normal computer.

## Install CIRI-long using pip

```bash
pip install CIRI-long
```
