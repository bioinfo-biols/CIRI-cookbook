# About

Manual of CIRI-full v2.0

If you have any questions, please contact 

- Yi Zheng @ Beijing Institutes of Life Science, Chinese Academy of Sciences.
- Email: zhengyi12@mails.ucas.ac.cn

CIRI-full is an accurate, high-throughput approach that uses both BSJ and reverse overlap (RO) features to reconstruct and quantify full-length circular RNAs from RNA-seq data sets. In CIRI-full, the BSJ feature is employed to detect cirexons and to determine the boundaries of circRNAs. The RO feature, deduced from the overlapped sequence of paired-end reads, is used to explore the detailed landscape within boundary sites. The alignments of both BSJ & RO merged reads will be visualized. The relative abundance of isoforms within one circRNA will be estimated according to the coverage and spliced events of BSJ & RO merged reads. 

# Installation

CIRI-full is developed in JAVA, and it can be performed in any system which has Java SE Runtime Environment.It requires:

```text
bwa:		A read mapping tool, which generates SAM file for CIRI-full, CIRI & CIRI-AS https://sourceforge.net/projects/bio-bwa/files/
CIRI2:		A circRNA detection tool  https://sourceforge.net/projects/ciri/
CIRI-AS:	A tool to detect cirexon and alternative splicing events in circRNAs https://sourceforge.net/projects/ciri/
```

CIRI2 and CIRI-AS are already packed with the CIRI-full software.

After downloading the CIRI-full package, you can extract it by typing:

```bash
unzip CIRI-full.zip	cd CIRI-full2. Preparation for running CIRI-full
```

Before running CIRI-full, you need to run CIRI and CIRI-AS to detect circRNAs and their associated BSJs and cirexons from your sequence data. 

Here is a recommend protocol to run CIRI and CIRI-AS: 

```bash
# Index the reference genome:
bwa index -a bwtsw reference.fa
# Split mapping using bwa-mem:
bwa mem -T 19 -t number_thread reference.fa read_1.fq read_2.fq > read.sam 
#  2.3 Running CIRI & CIRI-AS 
perl CIRI.pl -I read.sam -O prefix.ciri -F reference.fa -A annotation.gtf -T number_thread
perl CIRI_AS.pl -S read.sam -C prefix.ciri -F reference.fa -A annotation.gtf -O prefix -D yes
```

For detailed instructions on above tools, please read the manuals of bwa, CIRI and CIRI-AS.

# Running CIRI-full pipeline 

The CIRI-full Pipeline module is an automatic pipeline for detecting and reconstructing circRNAs. This pipeline includes `CIRI`, `CIRI-AS` and `CIRI-full` tools, which will finally generate reconstructed full-length circRNA sequences and the annotation of all identified circRNAs.

Before running the Pipeline module, please make sure that `bwa` is added to `$PATH`

The Pipeline module runs from a command line as follows:

```bash
java –jar CIRI-full.jar Pipeline [options]
```

Options:

```text
-1	reads1 of paired-end reads (required, equal length, fastq or fastq.gz format)
-2	reads2 of paired-end reads (required, equal length, fastq or fastq.gz format)
-r	reference genome in fasta format, the same file used in preparation step when building bwa index (required). 
-a	annotation file of reference genome in GTF format (optional).
-o	prefix of output files (optional, default: out)
-d	directory of output files (required)
-t	number of threads used in CIRI and bwa mem (optional, default: 1)	
-0	output all circRNAs including those with only one BSJ read support (optional, option for CIRI)
```

Four folders will be created under the dictionary set by `-d` option, `CIRI_output/`, `CIRI-AS_output/`, `CIRI-full_output/` and `sam/`, which contain the output files of `CIRI`, `CIRI-AS`, `CIRI-full` and `bwa`. 

For detailed information of these files, please refer to the following instructions.

# Running CIRI-full step-by-step

CIRI-full includes three modules, `RO1`, `RO2` and `Merge`. These modules should be performed sequentially in the following order: `RO1`, `RO2` and `Merge`.

## The RO1 module

This module is designed to identify 5’-RO feature on paired-end reads from RNA-seq data set and then, merge these RO containing paired-end reads into long single-end reads.

The RO1 module runs from a command line as follows:

```bash
java -jar CIRI_full.jar RO1 [options]
```

Options:

```text
-1	read1 of paired-end reads (required, equal length)	-2	read2 of paired-end reads (required, equal length)
-o	prefix of output files (optional,default: out)
-minM	sets the number of minimum 5’-RO length (optional, integer, default 13)
-minI	sets the minimum identity percentage of 5’-RO alignment (optional, default 95)
```

RO1 module will generate two output files:

```text
prefix_ro1_align.txt		
prefix_ro1.fq
```
		
Description of `prefix_ro1_align.txt`:

Each column gives the alignment information of each read pair which contain 5’-RO feature. 

- #read_id 
- #alignment_identity 
- #start_position_on_read1 
- #end_position_on_read1 
- #start_position_on_read2 
- #end_position_on_read2 
- #read_length

Description of `prefix_ro1.fq`

Read pairs with 5’-RO feature are merged into long sequences in FASTQ format. These sequences are taken as candidate RO merged-reads and will be filtered in the following steps.

## The RO2 module

The RO2 module is to analyze the alignment results of candidate RO merged-reads and screen out authentic ones for reconstructing full-length circRNAs.

Data preparation before running the RO2 module:

RO2 module filters RO merged-reads based on the SAM file generated by `bwa-mem`.

A recommended protocol for running bwa-mem:

```bash
bwa index -a bwtsw reference.fa	bwa mem -T 19 reference.fa prefix_ro1.fq > prefix_ro1.sam 
```

Note that `prefix_ro1.fq` file is the output file in the previous step (the RO1 module).

The RO2 module runs from a command line as follows:

```bash
java -jar CIRI_full.jar RO2 [options]Options:
```

Options:

```text
-r	reference genome in fasta format, the same file used in the preparation step when building bwa index (required). 
-s	SAM alignment of prefix.ro1.fq generated by bwa mem (required).
-l	the read length of given RNA-seq paired end data (required).
-range	maximum spanning distance of circRNAs on the reference(optional, integer, default 100000). 
-o	prefix of output files (required)
```

RO2 module will generate following output files:

```text
prefix_ro2.sam
prefix_ro2_info.list
```

Description of `prefix_ro2.sam`:

This file is the SAM alignment of authentic RO reads. 

Description of `prefix_ro2_info.list`:

This file gives the detailed alignment information of authentic RO reads. 

Columns are separated by tabs:

- #Read_ID
- #Chr
- #BSJ_position
- #Strand
- #Reconstructed_state
- #Cirexon
- #Mapping_order
- #Splice_site_state+ 
- #Splice_site_state-

#Splice_site_state+/- represents the mapping boundary deviation from the GT/AG splicing site, where -1 indicates that GT/AG splicing site cannot be detected on the current strand; positive value represents the distance between GT/AG splice site and split mapping position.

## The Merge module

The Merge module combines the results of RO2 and CIRI-AS to reconstruct full-length circRNAs. 

The Merge module runs from a command line as follows:

```bash
java –jar CIRI_full.jar Merge [options]
```

Options:

```text
-a	annotation file of the reference genome in GTF format (optional).
-c	output file of CIRI (required)
-as	output_all file generated in CIRI-AS (using -D yes argument. This file has a suffix “_jav.list” ) (required)
-ro	RO read information file (prefix_ro2_info.list) generated by RO2 module  (required)
-o	prefix of output files (required)
-r	reference genome file (in FASTA format) (required)
```

The Merge module will generate three output files.

```text
prefix_merge_circRNA_detail.anno  
```

Description of `prefix_merge_circRNA_detail.anno`

This file contains mapping information of BSJ reads (detected by CIRI) and RO merged-reads (detected by RO). Reads are clustered according to the BSJ position. Columns are separated by tabs:

- #BSJ
- #Chr
- #Start
- #End
- #GTF-annotated_exon
- #Cirexon
- #Coveage
- #BSJ_reads_information
- #RO_reads_information
- #Original_gene

# Running CIRI-vis

CIRI-vis is a tool for visualizing alignments of BSJ & RO merged reads and estimating the related abundance of isoforms according to the output of CIRI-full (prefix_merge_circRNA_detail.anno) or CIRI-AS (prefix_jav.list).

CIRI-vis.jar runs from a command line as follows:

```bash
java -jar CIRI-vis.jar [Options]
```

Options:

```text
-i		The path of input file of CIRI-vis. (required)
-l		The path of library length file. (required for isoform quantification)
-r		The path of reference genome sequence in FASTA format. (required for output circRNA sequence)
-list		The list of chosen circRNA BSJ. (optional)
-d		The dictionary of output. Default currentdir/stdir
-max	The maximum expression (BSJ reads number) of circRNA that displayed by CIRI-vis. Default 999999999
-min	The minimum expression (BSJ reads number) of circRNA that displayed by CIRI-vis. Default 10. **Note: please only use one of -min, -exp, -rank**
-rank	Only display the expression top X% of circRNA
-exp		Only display the top expression circRNA that contain X% of BSJ reads.
-iso		The maximum number of considering isoform, default 10. High value will make the quantification slower
```

CIRI-vis will output a set of pdf file, a “.list” file and a “.fa” file(if reference genome file is available) in a new created folder(set by “-d” parameter):

- One pdf file display circRNA isoforms on one BSJ.
- “.list” file shows detail information of each isoform. 
- “.fa” file shows the sequences of fully reconstructed circRNA isoforms.

Description of `prefix.list`:

This file gives the detailed information of circRNA isoforms. Columns are separated by tabs:

Columns|Description
---|---
1 | The name of pdf file.
2 | ID of the BSJ position of circRNA isoform in the pattern of "chr:start|end";
3 | chromosome of a predicted circRNA isoform
4 | start loci of a predicted circRNA isoform on the chromosome
5 | end loci of a predicted circRNA isoform on the chromosome
6 | circular junction read (also called as back-spliced junction read) count of a predicted circRNA
7 | the serial number of isoform in circRNA
8 | the estimate BSJ read count of this predicted isoform.
9 | the minimum length of this predicted isoform.
10 | whether this predicted isoform is fully reconstructed.
11 | The cirexon position in this predicted isoform, “0-0” represent for the breakpoint during reconstruction. 

Description of `prefix.fa`:

This FASTA format file will be generated if reference genome sequence is available. It  contains the sequence of fully reconstructed isoform. They were named in this format:

```text
>(Image_name)#(BSJ) length=(isoform_length) (isoform_BSJ_read_count)/(circRNA_BSJ_read_count)
```

If you want to display only a subset of circRNA, please use parameter “-list” to give CIRI-vis a list of BSJ position. The format should be like:

```text
chr10:74474869|74475660
chr8:141856359|141900868
```

**Notes:**
- IF you ran CIRI-full Pipeline in previous step, the input file will be named prefix_merge_circRNA_detail.anno under CIRI-full_output folder.
- IF you only ran CIRI-AS with '-d yes' parameter in previous step, the input file will be named prefix_jav.list under your CIRI-AS output folder.
- Library length file is necessary for isoform expression estimation. library length file will be prefix _library_length.list under your CIRI-AS output folder

# How to run the test data set using CIRI-full

Test data sets (FASTQ file, annotation file and reference sequence) are packaged with the CIRI-full software, which can be found in the “CIRI-full_test/“ folder. Temporary and final results are given in the “CIRI-full/test_output/” folder.

Here are the commands for running the test data sets: 

```bash
cd CIRI-full_v2.0/CIRI-full_test/
bwa index test_ref.fajava -jar ../CIRI-full.jar Pipeline -1 test_1.fq.gz -2 test_2.fq.gz -a test_anno.gtf -r test_ref.fa -d test_output/ -o testunset DISPLAY
java -jar ../CIRI-vis.jar -i test_output/CIRI-full_output/test_merge_circRNA_detail.anno -l ../CIRI-vis_test/test_library_length.list -r test_ref.fa –d test_output/CIRI-vis_out -min 1
```

If you want to run CIRI-full step by step, you can use the following commands:

```bash
cd CIRI-full_v2.0/CIRI-full_test/
mkdir test_output
bwa index test_ref.fa
bwa mem -T 19 test_ref.fa test_1.fq.gz test_2.fq.gz > test_output/test.sam 
perl ../bin/CIRI2.pl -I test_output/test.sam -O test_output/test.ciri -F test_ref.fa -A test_anno.gtf 
perl ../bin/CIRI_AS_v1.2.pl -S test_output/test.sam -C test_output/test.ciri -F test_ref.fa -A test_anno.gtf -O test_output/test -D yes
java -jar ../CIRI-full.jar RO1 -1 test_1.fq.gz -2 test_2.fq.gz -o test_output/test 
bwa mem -T 19 test_ref.fa test_output/test_ro1.fq > test_output/test_ro1.sam
java -jar ../CIRI-full.jar RO2 -r test_ref.fa -s test_output/test_ro1.sam -l 250 -o test_output/test
java -jar ../CIRI-full.jar Merge -c test_output/test.ciri -as test_output/test_jav.list -ro test_output/test_ro2_info.list -a test_anno.gtf -r test_ref.fa -o test_output/test
unset DISPLAY
java -jar ../CIRI-vis.jar -i test_output/CIRI-full_output/test_merge_circRNA_detail.anno -l ../CIRI-vis_test/test_library_length.list -r test_ref.fa -min 1
```

**Note: Please make sure you are using exactly the same version of genomic sequences and their annotations when running CIRI-full.**
