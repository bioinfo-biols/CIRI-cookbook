# Usage 3: Differential expression analysis

## Study without biological replicate

For sample without replicate, the differential expression & differential splicing analysis is 
performed using `CIRI_DE`

```
Usage:
  CIRI_DE [options] -n <control> -c <case> -o <out>

  <control>         CIRIquant result of control sample
  <case>            CIRIquant result of treatment cases
  <out>             Output file

Options (defaults in parentheses):

  -p                p value threshold for DE and DS score calculation (default: 0.05)
  -t                numer of threads (default: 4)

Example usage:
  CIRI_DE -n control.gtf -c case.gtf -o CIRI_DE.tsv
```

The output format `CIRI_DE` is in the format below:

| column | name | description |
|--------|------|-------------|
| 1 | circRNA_ID | circRNA identifier |
| 2 | Case_BSJ | number of BSJ reads in case |
| 3 | Case_FSJ | number of FSJ reads in case |
| 4 | Case_Ratio | junction ratio in case |
| 5 | Ctrl_BSJ | number of BSJ reads in control |
| 6 | Ctrl_FSJ | number of FSJ reads in control |
| 7 | Ctrl_Ratio | junction ratio  in control |
| 8 | DE_score | differential expression score |
| 9 | DS_score | differential splicing score |

## Study with biological replicates

For study with biological replicates, a customed analysis pipeline of edgeR is recommended and 
we provide `prep_CIRIquant` to generate matrix of circRNA expression level / junction ratio and `CIRI_DE_replicate` 
for DE analysis

**Step1**: Prepare CIRIquant output files

One should provide a text file listing sample information and path to CIRIquant output GTF files

```
CONTROL1 ./c1/c1.gtf C 1
CONTROL2 ./c2/c2.gtf C 2
CONTROL3 ./c3/c3.gtf C 3
CASE1 ./t1/t1.gtf T 1
CASE2 ./t2/t2.gtf T 2
CASE3 ./t3/t3.gtf T 3
```

The first three columns is required by default. For paired samples, you could also add a column of subject name.

| column | description |
|--------|-------------|
| 1 | sample name |
| 2 | path to CIRIquant output gtf |
| 3 | group ("C" for control, "T" for treatment) |
| 4 | subject (optional, only for paired samples) |

**Note: If you are planning to use CIRI_DE for differential expression, then group name in column 3 must be either "C" or "T".**

Then, run `prep_CIRIquant` to summarize the circRNA expression profile in all samples

```
Usage:
  prep_CIRIquant [options]

  -i                the file of sample list
  --lib             where to output library information
  --circ            where to output circRNA annotation information
  --bsj             where to output the circRNA expression matrix
  --ratio           where to output the circRNA junction ratio matrix

Example:
  prep_CIRIquant -i sample.lst \
                 --lib library_info.csv \
                 --circ circRNA_info.csv \
                 --bsj circRNA_bsj.csv \
                 --ratio circRNA_ratio.csv
```

These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR 
(using the DESeqDataSetFromMatrix and DGEList functions, respectively).

**Step2**: Prepare StringTie output

The output of StringTie should locate under `output_dir/gene/prefix_out.gtf`. You need to use 
[prepDE.py](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py) from stringTie to
generate the gene count matrix for normalization.

For example, one can provide a text file `sample_gene.lst` containing sample IDs and path to StringTie outputs:

```text
CONTROL1 ./c1/gene/c1_out.gtf
CONTROL2 ./c2/gene/c2_out.gtf
CONTROL3 ./c3/gene/c3_out.gtf
CASE1 ./t1/gene/t1_out.gtf
CASE2 ./t2/gene/t2_out.gtf
CASE3 ./t3/gene/t3_out.gtf
```

Then, run `prepDE.py -i sample_gene.lst` and use `gene_count_matrix.csv` generated under current working directory 
for further analysis.

**Step3**: Differential expression analysis

For differential analysis using `CIRI_DE_replicate`, you need to install a R environment and `edgeR` package from Bioconductor.

```bash
usage: CIRIquant_DE_replicate [-h] --lib FILE --bsj FILE --gene FILE --out
                              FILE --out2 FILE

optional arguments:
  -h, --help   show this help message and exit
  --lib FILE   library information
  --bsj FILE   circRNA expression matrix
  --gene FILE  gene expression matrix
  --out FILE   output result of circRNA differential expression analysis
  --out2 FILE  output result of gene differential expression analysis

Example:
  CIRI_DE_replicate \
          --lib  library_info.csv \
          --bsj  circRNA_bsj.csv \
          --gene gene_count_matrix.csv \
          --out  circRNA_de.tsv \
          --out2 gene_de.tsv
```

Please be noted that the output results is **unfiltered**,
and you could apply a more stringent filter on expression values to get a more convincing result.
