# Usage

## Basic usage

```
usage: CIRI-long [-h] [-v] {call,collapse} ...

positional arguments:
  {call,collapse}  commands

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show program's version number and exit
```

CIRI-long have two main functions, including (1) candidate circRNAs identification and (2) isoform collapsing. 

## Step1. circRNA identification

### Basic options

```
usage: CIRI-long call [-h] [-i READS] [-o DIR] [-r REF] [-p PREFIX] [-a GTF] [--canonical] [-t INT] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -i READS, --in READS  Input reads.fq.gz
  -o DIR, --out DIR     Output directory, default: ./
  -r REF, --ref REF     Reference genome FASTA file
  -p PREFIX, --prefix PREFIX
                        Output sample prefix, (default: CIRI-long)
  -a GTF, --anno GTF    Genome reference gtf, (optional)
  -c CIRC, --circ CIRC  Additional circRNA annotation in bed/gtf format,
                        (optional)
  --canonical           Use canonical splice signal (GT/AG) only, default: True)
  -t INT, --threads INT
                        Number of threads, (default: use all cores)
  --debug               Run in debugging mode, (default: False)
```

**NOTE**:
- A bwa index for reference genome is required, please use `bwa index` command to generate bwa index before running CIRI-long.

### Example Usage:

Demo dataset can be downloaded from the [GitHub release](https://github.com/bioinfo-biols/CIRI-long/releases)

```
# Download demo dataset
wget https://github.com/bioinfo-biols/CIRI-long/releases/download/v0.6-alpha/CIRI-long_test_data.tar.gz

# Decompress demo dataset
tar zxvf CIRI-long_test_data.tar.gz
cd test_data

# Build bwa index before running CIRI-long
bwa index -a bwtsw mm10_chr12.fa mm10_chr12.fa

# Run CIRI-long to identify circular reads from sequencing reads
CIRI-long call -i test_reads.fa \
               -o ./test_call \
               -r mm10_chr12.fa \
               -p test \
               -a mm10_chr12.gtf \
               -t 8
```

### Output Files

The output directory should have the following structure:

```bash
test_call
├── test.cand_circ.fa
├── test.json
├── test.log
├── test.low_confidence.fa
└── tmp
    ├── ss.idx
    ├── test.ccs.fa
    └── test.raw.fa

1 directory, 7 files
```

## Step2. isoform collapse

### Basic Options

```
usage: CIRI-long collapse [-h] [-i LIST] [-o DIR] [-p PREFIX] [-r REF] [-a GTF] [--canonical] [-t INT] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -i LIST, --in LIST    Input list of CIRI-long results
  -o DIR, --out DIR     Output directory, default: ./
  -p PREFIX, --prefix PREFIX
                        Output sample prefix, (default: CIRI-long)
  -r REF, --ref REF     Reference genome FASTA file
  -a GTF, --anno GTF    Genome reference gtf, (optional)
  -c CIRC, --circ CIRC  Additional circRNA annotation in bed/gtf format,
                        (optional)
  --canonical           Use canonical splice signal (GT/AG) only, default: True)
  -t INT, --threads INT
                        Number of threads, (default: use all cores)
  --debug               Run in debugging mode, (default: False)
```

One should provide a text file listing sample name and path to CIRI-long output files `*.cand_circ.fa`, seperated by space.

```text
sample1_name /path/to/sample1/cand_circ.fa
sample2_name /path/to/sample2/cand_circ.fa
```

### Example Usage

For exmaple, you can create a file name `test.lst` with the following content:

```text
test ./test_call/test.cand_circ.fa
```

Then run `CIRI-long collapse` to aggregate results from one or multiple samples.

```bash
 CIRI-long collapse -i ./test.lst \
                    -o ./test_collpase \
                    -p test \
                    -r ./mm10_chr12.fa \
                    -a ./mm10_chr12.gtf \
                    -t 8
```

### Output Files

The output directory should have the following structure:

```
test_collpase
├── test_collpase.expression
├── test_collpase.info
├── test_collpase.log
├── test_collpase.reads
└── tmp
    ├── ss.idx
    └── test_collpase.corrected.pkl

1 directory, 6 files
```

### Output Format

**The main output**

The main output of CIRI-long is a GTF file (e.g. `test_collpase.info`), that contains detailed information of circRNAs and annotation of circRNA back-spliced regions in the attribute columns

Description of each columns's value

| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | source | CIRI-long |
| 3 | type | circRNA |
| 4 | start | 5' back-spliced junction site |
| 5 | end | 3' back-spliced junction site |
| 6 | score | Number of total supported reads |
| 7 | strand | strand information |
| 8 | . | . |
| 9 | attributes | attributes seperated by semicolon |

The attributes containing several pre-defined keys and values:

| key| description|
|----|------------|
|circ_id | name of circRNA|
|splice_site | splicing signal of candidate circRNAs and numbers indicating shifted bases of aligned and annotated splice site. (e.g. AG-GT\|0-5) |
|equivalent_seq | equivalent sequence of splice site |
|circ_type| circRNA types: exon/intron/intergenic|
|circ_len| length of the major isoform of circRNA|
|isoform| structure of isoforms, isoforms are seperated by "\|" and circular exons are seperated by "," (e.g. 11627815-111627914,111628190-111628302\|11627815-111628302) |
|gene_id| ensemble id of host gene |
|gene_name |HGNC symbol of host gene|
|gene_type | type of host gene in the annotation gtf file|

**Expression matrix**

`test_collpase.expression` contains the summarized expression level of circRNAs in all samples in `tsv` format.

## Using non-canonical splice signals

If you would like to use other splice signals, please modify the dict `SPLICE_SIGNAL` in [align.py](https://github.com/Kevinzjy/CIRI-long/blob/master/CIRI/align.py#L34) in format: {(5'SS, 3'SS): Priority}

Default configuration:

```Python
SPLICE_SIGNAL = {
    ('GT', 'AG'): 0,  # U2-type
    ('GC', 'AG'): 1,  # U2-type
    ('AT', 'AC'): 2,  # U12-type
    ('GT', 'AC'): 2,  # U12-type
    ('AT', 'AG'): 2,  # U12-type
}
```

## Using additional circRNA annotations

From version `v1.0.2`, CIRI-long call also provide additional circRNA annotations in BED/GTF format for BSJ correction with `--circ` option. CircRNA annotations can be downloaded from [circAtlas](http://circatlas.biols.ac.cn/) or other databases. The GTF-format output of CIRIquant is also supported.

**NOTE: If using results from other tools/databases, please make sure the coordinate system is compatible with our CIRI-series tools**:

> The coordinate system of circRNAs is different in most circRNA tools. For instance, if a circRNA is derived from chr1:1000-2000, it should be reported as chr1:1000-2000 in CIRI-series and some tools (DCC/KNIFE/Mapsplice), but reported as chr1:999-2000 in other tools (CIRCexplorer2/UROBORUS/circRNA_finder/find_circ).
> 
> Thus, if you want to use circRNAs identified from tools in the latter group, you need to add 1 extra base to the start coordinate of circRNAs (the position with smaller coordinate regardless of the strand information), then use the altered coordinates as input.
