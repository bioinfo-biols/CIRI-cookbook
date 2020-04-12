# Usage 1: circRNA quantifcation

## Basic options

```text
Usage:
  CIRIquant [options] --config <config> -1 <m1> -2 <m2>

  <config>          Config file
  <m1>              Input mate1 reads (for paired-end data)
  <m2>              Input mate2 reads (for paired-end data)


Options (defaults in parentheses):

  -v                Run in verbose mode
  -o, --out          Output directory (default: current directory)
  -e, --log         Specific log file (default: sample_prefix.log)
  -p, --prefix      Output sample prefix (default: input sample name)
  -t, --threads     Number of CPU threads to use (defualt: 4)

  --bed             User provided Back-Spliced Junction Site in BED format
  --circ            User provided circRNA prediction results
  --tool            User provided tool name for circRNA prediction (Only CIRI2 is supported)

  --RNaseR          CIRIquant output file of RNase R data (required for RNase R correction)
  --bam             Specific hisat2 alignment bam file against reference genome
  --no-gene         Skip StringTie estimation of gene abundance
```

**NOTE**: 
- For now, --circ and --tool options can only parse CIRI2 results.
- Gene expression values are needed for normalization, do not use `--no-gene` if you need to run DE analysis afterwards. 

## Example config

A YAML-formated config file is needed for CIRIquant to find software and reference needed.
A valid example of config file is demonstrated below.

```YAML
// Example of config file
name: hg19
tools:
  bwa: /home/zhangjy/bin/bwa
  hisat2: /home/zhangjy/bin/hisat2
  stringtie: /home/zhangjy/bin/stringtie
  samtools: /home/zhangjy/bin/samtools

reference:
  fasta: /home/zhangjy/Data/database/hg19.fa
  gtf: /home/zhangjy/Data/database/gencode.v19.annotation.gtf
  bwa_index: /home/zhangjy/Data/database/hg19/_BWAtmp/hg19
  hisat_index: /home/zhangjy/Data/database/hg19/_HISATtmp/hg19
```

Key | Description
----|-------------
name| the name of config file
bwa | the path of `bwa`
hisat2 | the path of `hisat2`
stringtie | the path of `stringite`
samtools | the path of `samtools`, samtools version below 1.3.1 is not supported
fasta | reference genome fasta, a fai index by `samtools faidx` is also needed under the same directory
gtf | annotation file of reference genome
bwa_index | prefix of BWA index for reference genome
hisat_index | prefix of HISAT2 index for reference genome

For quantification of user-provided circRNAs, a list of junction sites in bed format is required, for example:

```
chr1    10000   10099   chr1:10000|10099    .   +
chr1    31000   31200   chr1:31000|31200    .   -
```

## Example Usage

Test data set can be retrived from [test_data.tar.gz](https://github.com/Kevinzjy/CIRIquant/releases/download/v0.2.0/test_data.tar.gz), you can replace the path of required software in
 the `chr1.yml` with your own version

```
tar zxvf test_data.tar.gz
cd test_data/quant
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          --no-gene \
          -o ./test \
          -p test
```

The output file `test.gtf` should be located under `test_data/quant/test`  

The demo dataset should take approximately 5 minutes on a personal computer. It has been tested on 
my PC with Intel i7-8700 processor and 16G of memory, running Ubuntu 18.04 LTS.

## Output format

The main output of CIRIquant is a GTF file, that contains detailed information of 
BSJ and FSJ reads of circRNAs and annotation of circRNA back-spliced regions in the attribute columns

Description of each columns's value

| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | source | CIRIquant |
| 3 | type | circRNA
| 4 | start | 5' back-spliced junction site |
| 5 | end | 3' back-spliced junction site |
| 6 | score | CPM of circRNAs (#BSJ / #Mapped reads)
| 7 | strand | strand information
| 8 | . | .
| 9 | attributes | attributes seperated by semicolon |

The attributes containing several pre-defined keys and values:

| key | description|
| --- | -----------|
| circ_id | name of circRNA |
| circ_type | circRNA types: exon / intron / intergenic |
| bsj | number of bsj reads |
| fsj | number of fsj reads |
| junc_ratio | circular to linear ratio: 2 * bsj / ( 2 * bsj + fsj)
| rnaser_bsj | number of bsj reads in RNase R data (only when --RNaseR is specificed)
| rnaser_fsj | number of fsj reads in RNase R data (only when --RNaseR is specificed)
| gene_id | ensemble id of host gene |
| gene_name | HGNC symbol of host gene |
| gene_type | type of host gene in gtf file |
