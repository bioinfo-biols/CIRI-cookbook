# About

CIRI-vis is a tool for visualizing alignments of BSJ & RO merged reads and estimating the related abundance of isoforms according to the output of CIRI-full (prefix_merge_circRNA_detail.anno) or CIRI-AS (prefix_jav.list).

![CIRI-vis.png](https://github.com/bioinfo-biols/Zhaolab/blob/master/CIRI-vis.png?raw=true)

## Author

Authors: Yi Zheng(zhengyi@biols.ac.cn), Fangqing Zhao(zhfq@biols.ac.cn)

Maintainer: Yi Zheng

## Release Notes

- Version 1.4

# Installation

## Prerequisites

```
Softwares:
	JavaSE >= 1.6
    bwa
	CIRI2
    CIRI-AS
    CIRI-Full
```

## Install CIRI-vis

CIRI-vis is developed in JAVA, and it can be performed in any system which has Java SE Runtime Environment.

CIRI-vis is already packed with the CIRI-full software under /bin.
After downloading the CIRI-full package, you can extract it by typing:

```bash
unzip CIRI-full_v2.0.zip
```

you can find CIRI-vis.jar in the folder.

# Commands and Arguments

## Running CIRI-vis

Input file requirements:

IF you runned CIRI-full Pipeline in previous step, the input file will be named:XXX_merge_circRNA_detail.anno under CIRI-full_output folder 

IF you only run CIRI-AS with '-d yes' parameter in previous step, the input file will be named XXX_jav.list under your CIRI-AS output folder 

Library length file is nessaracy for isoform expression estimation. library length file will be XXX_library_length.list under your CIRI-AS output folder

CIRI-vis.jar runs from a command line as follows:

```text
Usage: java -jar CIRI-vis.jar [Options]

Options:

	-i			The path of input file of CIRI-vis. (required)
	-l			The path of library length file. (required for isoform quantification)
	-r			The path of reference genome sequence in FASTA format. (required for output circRNA sequence)
	-list		The list of choosen circRNA BSJ.(It is needed when more than one sample)
	-d			The dictionary of output. Default currentdir/stdir
	-o			The prefix of output. Default stout (optional)
	-type		The format of figure. you can select pdf or svg or both. Default pdf (optional)
	-max		The maximum expression (BSJ reads number) of circRNA that displayed by CIRI-vis. Default 999999999 (optional)
	-min		The minimum expression (BSJ reads number) of circRNA that displayed by CIRI-vis. Default 5. Note: please only use one of -min, -exp, -rank (optional)
	-rank		Only display the expression top X% of circRNA (optional)
	-exp		Only display the top expression circRNA that contain X% of BSJ reads. (optional)
	-iso		The maximum number of considering isoform, default 10. High value will make the quantification slower (optional)
	-ran		Set random seed, default 0.(optional)
```

Examples:

For one sample

```bash
java -jar CIRI-vis.jar -i A_merge_circRNA_detail.anno -l A_library_length.list -r Ref.fa -d out -o prefix
```

For more than one samples

```bash
java -jar CIRI-vis.jar -i A_merge_circRNA_detail.anno B_merge_circRNA_detail.anno -l A_library_length.list B_library_length.list -r Ref.fa -d out -o prefix -list test.txt
```

The format of the "test.txt" (input of -list) should be:

```text
chr10:74474869|74475660
chr8:141856359|141900868	
...
```

## Output files

Description of `prefix.list`:

This file gives the detailed information of circRNA isoforms. Columns are separated by tabs:

```text
Image_ID:		The name of pdf file.
Circle_ID:		ID of the BSJ position of circRNA isoform in the pattern of "chr:start|end";
Chr: 			chromosome of a predicted circRNA isoform
start:			start loci of a predicted circRNA isoform on the chromosome
end:			end loci of a predicted circRNA isoform on the chromosome
total_exp:		circular junction read (also called as back-spliced junction read) count of a predicted circRNA
isoform_number: the serial number of isoform in circRNA
isoform_exp:	the estimate BSJ read count of this predicted isoform.
isoform_length: the minimum length of this predicted isoform.
isoform_state:	whether this predicted isoform is fully reconstructed.
strain: 		strain of circRNA (+/-)
gene_id:		gene name that the circRNA located in
isoform_cirexon: The cirexon position in this predicted isoform, "0-0" represent for the breakpoint during reconstruction. 
```

When more than one sample

```text
Image_ID:		The name of pdf file.
Sample_name:	 Name of input sample (appear only when more than one sample)
Circle_ID:		ID of the BSJ position of circRNA isoform in the pattern of "chr:start|end";	
Chr: 			chromosome of a predicted circRNA isoform
start:			start loci of a predicted circRNA isoform on the chromosome
end:			end loci of a predicted circRNA isoform on the chromosome
total_exp:		circular junction read (also called as back-spliced junction read) count of a predicted circRNA
isoform_number: the serial number of isoform in circRNA
isoform_exp:	the estimate BSJ read count of this predicted isoform.
isoform_length: the minimum length of this predicted isoform.
estimated_isoform_read_count:	The estimated total number of read that on this isoform (including BSJ and nonBSJ reads),
isoform_state:	whether this predicted isoform is fully reconstructed.
strain: 		strain of circRNA (+/-)
gene_id:		gene name that the circRNA located in
isoform_cirexon: The cirexon position in this predicted isoform, "0-0" represent for the breakpoint during reconstruction. 
```

Description of `prefix.fa`:

This FASTA format file will be generated if reference genome sequence is available. It  contains the sequence of fully reconstructed isoform. They were named in this format:

```text
>(Image_name)#(BSJ) length=(isoform_length) (isoform_BSJ_read_count)/(circRNA_BSJ_read_count)
```

# Example Usage

Test data sets (FASTQ file, annotation file and reference sequence) are packaged with the CIRI-full software, which can be found in the CIRI-full_test/ folder. Temporary and final results are given in the CIRI-full/test_output/ folder.

Here are the commands for running the test data sets: 

```bash
cd CIRI-full_v2.0/CIRI-full_test/
bwa index test_ref.fa
java -jar ../CIRI-full.jar Pipeline -1 test_1.fq.gz -2 test_2.fq.gz -a test_anno.gtf -r test_ref.fa -d test_output/ -o test
unset DISPLAY
java -jar CIRI-vis.jar -i test_output/CIRI-full_output/test_merge_circRNA_detail.anno -l ../CIRI-vis_test/test_library_length.list -r test_ref.fa -d test_output/CIRI-vis_out -min 1
```
