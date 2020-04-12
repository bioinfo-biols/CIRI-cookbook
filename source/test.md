# Test data

A testing dataset can be downloaded from [Github](https://github.com/Kevinzjy/CIRIquant/releases/download/v0.1.0/test_data.tar.gz)

Folder `quant` contain the test dataset for circRNA quantification

## 1. Generate hisat2 and bwa index

```bash
bwa index -a bwtsw -p chr1.fa chr1.fa
hisat2-build ./chr1.fa ./chr1.fa
```

## 2. Edit the configuration in chr1.yml

Change the path of bwa/hisat2/stringtie/samtools to your own path

## 3. Run test dataset

Test data set can be retrived under `test_data/quant` folder, you can replace the path of required software in the `chr1.yml` with your own version

```bash
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          --no-gene \
          -o ./test \
          -p test
```

The structure of output directory `./test` should be like this:

```text
test
├── align
│   ├── test.bam
│   ├── test.sorted.bam
│   └── test.sorted.bam.bai
├── circ
│   ├── test.ciri
│   ├── test.ciri.bed
│   ├── test_denovo.bam
│   ├── test_denovo.sorted.bam
│   ├── test_denovo.sorted.bam.bai
│   ├── test_index.1.ht2
│   ├── test_index.2.ht2
│   ├── test_index.3.ht2
│   ├── test_index.4.ht2
│   ├── test_index.5.ht2
│   ├── test_index.6.ht2
│   ├── test_index.7.ht2
│   ├── test_index.8.ht2
│   ├── test_index.fa
│   └── test_unmapped.sam
├── CIRIerror.log
├── test.bed
├── test.gtf
└── test.log
```

You can check the output in `./test/test.gtf` and compare it to `test.gtf` provided in test_data

## Differential expression analysis

Folder `DE` contain the test dataset for differential expression analysis

```bash
cd test_data/DE

# Test for DE-score and DS-score calculation
CIRI_DE -n ctrl.gtf \
        -c case.gtf \
        -o CIRI_DE.csv

# Test for RNase R correction
CIRI_DE -n ctrl_corrected.gtf \
        -c case_corrected.gtf \
        -o CIRI_DE_corrected.csv
```
