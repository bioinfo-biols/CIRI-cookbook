# Usage 1: circRNA quantifcation

## Basic options

```
usage: CIRIquant [-h] [--config FILE] [-1 MATE1] [-2 MATE2] [-o DIR]
                 [-p PREFIX] [-t INT] [-a INT] [-l INT] [--ciri3] [-v]
                 [--version] [-e LOG] [--bed FILE] [--circ FILE] [--tool TOOL]
                 [--RNaseR FILE] [--bam BAM] [--no-gene] [--no-fsj]
                 [--bsj-file FILE]

optional arguments:
  -h, --help            show this help message and exit
  --config FILE         Config file in YAML format
  -1 MATE1, --read1 MATE1
                        Input mate1 reads (for paired-end data)
  -2 MATE2, --read2 MATE2
                        Input mate2 reads (for paired-end data)
  -o DIR, --out DIR     Output directory, default: ./
  -p PREFIX, --prefix PREFIX
                        Output sample prefix, default: input sample name
  -t INT, --threads INT
                        Number of CPU threads, default: 4
  -a INT, --anchor INT  Minimum anchor length for junction alignment, default:
                        5
  -l INT, --library-type INT
                        Library type, 0: unstranded, 1: read1 match the sense
                        strand,2: read1 match the antisense strand, default: 0
  -v, --verbose         Run in debugging mode
  --version             show program's version number and exit
  -e LOG, --log LOG     Log file, default: out_dir/prefix.log
  --bed FILE            bed file for putative circRNAs (optional)
  --circ FILE           circRNA prediction results from other softwares
  --tool TOOL           circRNA prediction tool, required if --circ is
                        provided
  --RNaseR FILE         CIRIquant result of RNase R sample
  --bam BAM             hisat2 alignment to reference genome
  --no-gene             Skip stringtie estimation for gene abundance
  --no-fsj              Skip FSJ extraction to reduce run time
  --bsj-file FILE       output BSJ read IDs to file (optional)
```

**NOTE**:
- For now, --circ and --tool options support results from `CIRI2` / `CIRCexplorer2` / `DCC` / `KNIFE` / `MapSplice` / `UROBORUS` / `circRNA_finder` / `find_circ`
- For tools like `DCC` and `circRNA_finder`, please manually remove duplicated circRNAs with same junction postion but have opposite strands.
- Gene expression values are needed for normalization, do not use `--no-gene` if you need to run DE analysis afterwards.

## Example YAML config

A YAML-formated config file is needed for CIRIquant to find software and reference needed.

A valid example of minimal config file:

```YAML
reference:
  fasta: /home/zhangjy/Data/database/hg19.fa
  gtf: /home/zhangjy/Data/database/gencode.v19.annotation.gtf
  bwa_index: /home/zhangjy/Data/database/hg19/_BWAtmp/hg19
  hisat_index: /home/zhangjy/Data/database/hg19/_HISATtmp/hg19
```

An example of supported config file:

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
name| the name of config file (optional)
bwa | the path of `bwa` (optional, defaults to bwa in $PATH)
hisat2 | the path of `hisat2` (optional, defaults to hisat2 in $PATH)
stringtie | the path of `stringite` (optional, defaults to stringtie in $PATH)
samtools | the path of `samtools`, samtools version below 1.3.1 is not supported (optional, defaults to samtools in $PATH)
fasta | reference genome fasta, a fai index by `samtools faidx` is also needed under the same directory
gtf | annotation file of reference genome in GTF/GFF3 format
bwa_index | prefix of BWA index for reference genome
hisat_index | prefix of HISAT2 index for reference genome

## Example circRNA bed file

For quantification of user-provided circRNAs, a list of junction sites in bed format is required, the 4th column must be in "chrom:start|end" format. For example:

```text
chr1    10000   10099   chr1:10000|10099    .   +
chr1    31000   31200   chr1:31000|31200    .   -
```

## Example Usage

### Recommended: Predict circRNAs using CIRI2 (packaged in CIRIquant)

```bash
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          -o ./test \
          -p test
```

### Quantify circRNAs using provided BED format input

```
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          -o ./test \
          -p test \
          --bed your_circRNAs.bed
```

### Quantify circRNAs using results from other tools

For example, if you have `find_circ` results of predicted circRNAs.

```
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          -o ./test \
          -p test \
          --circ find_circ_results.txt \
          --tool find_circ
```

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
