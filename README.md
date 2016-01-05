# TIS_tools
Software tools for processing sequencing data from transposon insertion sequencing (TIS) experiments

## Introduction

Transposon insertion sequencing (TIS) harnesses the power and throughput of next generation sequencing to allow for genome-wide analysis of microbial gene function and pathways. Several methods for performing TIS library preparation and sequencing have been described (TnSeq, INSeq, HITS, TraDIS) and software tools such as INSeq_analysis and ESSENTIALS have been created to analyze sequencing results. TIS_tools is a set of perl scripts developed as a bridge between sequencing and analyisis by performing read alignments, filtering alignments, and outputting read count data in a format suitable for downstream analysis. 

## Scripts:

### INSeq_read_preprocess.pl 

For sequence reads generated using the "INSeq" transposon insertion sequencing protocol, performs the following processing steps:
  
  1. Separate multiplexed libraries by "barcode" sequence
  2. Identify and trim transposon sequence from reads
  3. Align reads to reference sequence(s)
  4. Assign aligments to positions in reference sequence
  5. Read count correction and filtering
  6. Output read counts
  
  **Prerequisites:**   
  [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (version 0.12.8 or higher)
  
  **Usage:**   
  `perl INSeq_read_preprocess.pl [options] -r <reads.fastq> -c <barcodes.txt> -g <reference_sequence.fasta>`

  **Required arguments:**  
  `-r` Read sequence file, in fastq format. Can remain gzipped.  
  `-c` Tab-delimited text file of barcodes in the following format, one barcode per line:  
  `pool_id<tab>barcode`  
  *Optional:* If you are re-processing already de-barcoded and transposon-trimmed reads or if reprocessing previously-performed bowtie alignments, you can also include paths to either the processed read output files with or without paths to the bowtie alignment files. To include the sorted reads file (as previously output by this script): `pool_id<tab>barcode<tab>/path/to/reads.fasta`. To also include the bowtie alignment (as previously output by this program): `pool_id<tab>barcode<tab>/path/to/reads.fasta<tab>/path/to/alignment.bowtie`.  
  `-g` Reference genome sequence file, in fasta format. For best results, this file must include sequences for ALL genetic material present in the organism (i.e. chromosomes and plasmids).
  
  **Encouraged, but optional:**  
  `-w` Sequence of parent strain that was used ot produce transposon mutagenesis library, in fasta format. OK if this is a multi-contig draft assembly.
  
  **Other options:**  
  `-t` Transposon sequence (default: "ACAGGTTG")  
  `-m` Number of transposon sequence mismatches allowed (default: 1)  
  `-o` Output format for count files. Possible options include:
  + 'wiggl': Files will be in wiggle format that can be used as input to ESSENTIALS (default)
  + 'inseq': Files will be in in the format of *_bowtiemap_processed.txt files produced by Goodman et al's INSeq_analysis software
  + 'essen': Files will be in the format of allta_split_*.counts.txt files produced by ESSENTIALS
  `-s` Minimum number of total reads (left flank + right flank) at an insertion site required for output (default: 3)  
  `-d` Minimum difference in read counts betwen sides to trigger investigation of site flanks. A higher number may result in faster processing, but lower sensitivity for identifying read count discrepancies. Minimum value is 1. (default: 1)  
  `-h` Ignore "shifted" reads, i.e. reads that align perfectly without mismatches, but one base upstream or downstream from a "TA" insertion site. (default: "shifted" reads will be added to a flank's total read count if addition of the reads results in the two flanks having a more balanced read count)  
  `-b` Path to directory containing bowtie and bowtie-build. For example: /Users/myname/applications/bowtie_folder (default: assumes this directory is in your PATH)  
  `-p` Number of parallel threads to run (default: 1)

**Output files:**    
