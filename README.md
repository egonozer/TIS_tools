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

  **Required Arguments:**
  `-r` Read sequence file, in fastq format. Can remain gzipped.
  `-c` Tab-delimited text file of barcodes in the following format, one barcode per line:
  `pool_id<tab>barcode`
  Optional: If you are re-processing already de-barcoded and transposon-trimmed reads or if reprocessing previously-performed bowtie alignments, you can also include paths to either the processed read output files with or without paths to the bowtie alignment files. To include the sorted reads file (as previously output by this script): `pool_id<tab>barcode<tab>/path/to/reads.fasta`. To also include the bowtie alignment (as previously output by this program): `pool_id<tab>barcode<tab>/path/to/reads.fasta<tab>/path/to/alignment.bowtie`.
  `-g` Reference genome sequence file, in fasta format. For best results, this file must include sequences for ALL genetic material present in the organism (i.e. chromosomes and plasmids).
  
