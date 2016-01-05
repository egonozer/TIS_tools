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

**Program output:**  
As the script runs, several statistics will be output to the screen. If you wish to capture these results, simply redirect STDOUT to a file like so:  
`perl INSeq_read_preprocess.pl -r reads.fastq -c barcodes.txt -g reference_sequence.fasta > stats.txt`  
or:  
`perl INSeq_read_preprocess.pl -r reads.fastq -c barcodes.txt -g reference_sequence.fasta | tee stats.txt`  

Barcode sorting, trimming, and transposon trimming results will take the following format:  
```
Total reads: xxx
Total reads with barcode: xxx
Total reads with transposon: xxx
Transposons with 0 mismatch(es): xxx
Transposons with 1 mismatch(es): xxx
pool	bc	#_w_bc	%_w_bc	#_w_tn %_w_tn
...
```
`Total reads` is the total number of reads found in the input read fastq file  
`Total reads with barcode` is the total number of reads found to have a recognizable barcode sequence at the start of the read sequence  
`Total reads with transposon` is the total number of reads (from those found to have a barcode sequence) that had the given transposon sequence at the expected position (16 or 17 bp from the 5' end after removal of barcode sequence)  
`pool`: Pool ID  
`bc`: Barcode  
`#_w_bc`: Number of reads with this barcode  
`%_w_bc`: Percentage of total reads with this barcode  
`#_w_tn`: Number of reads with this barcode found to have the given transposon sequence in the expected position  
`%_w_tn`: Percentage of reads with this barcode found to have the given transposon sequence in the expected position

**Output files:**    
Files beginning with "temp_" can be safely deleted.  
Files will start with a prefix corresponding to the pool_id and barcode given in the barcode text file above. For example, if one of the pool IDs was "Input_1" identified with barcode "TATA", then all files relating to this pool will be prefixed with "Input_1_TATA".  
- `<prefix>.reads.fasta`: Fasta formatted file of all sorted, de-barcoded, and transposon trimmed read sequences associated with this pool. 
- `<prefix>.bowtie`: Alignment file produced by bowtie
- `<prefix>.<reference_sequence>.wiggle`: Read alignment counts for each of the individual records that was present in the reference genome sequence file. For example, if the file given to the `-r` option contained one chromosomal sequence record and one plasmid sequence record there should be two of these files per barcode, each named according to the sequence record name and with coordinates corresponding to positions along the indicated sequence. If outputting in the default "wiggle" format, each line of the file will contain a sequence position (1-based) corresponding to a transposon insertion site and the number of reads aligning to that position separated by a tab. When a particular position is listed twice, the first listing is the number of reads aligning to the left flank (upstream) of the insertion site and the second listing is the number of reads aligning to the right flank (downstream) of the insertion site.

