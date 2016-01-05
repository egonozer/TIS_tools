# TIS_tools
Software tools for processing sequencing data from transposon insertion sequencing (TIS) experiments

## Introduction

Transposon insertion sequencing (TIS) harnesses the power and throughput of next generation sequencing to allow for genome-wide analysis of microbial gene function and pathways. Several methods for performing TIS library preparation and sequencing have been described (TnSeq, INSeq, HITS, TraDIS) and software tools such as INSeq_analysis and ESSENTIALS have been created to analyze sequencing results. TIS_tools were developed to bridge between sequencing and analyisis by performing read alignments, filtering alignments, and outputting read count data in a format suitable for downstream analysis.

## Concept / Methods

Analysis of transposon insertion sequencing experiments depends on accurate determination of the number of sequencing reads aliging to a particular transposon insertion location to approximate the numbers and identities of particular transposon mutants present before and after application of a selective condition. There are several potential reasons why aligned read counts could be inaccurate:

1. Sequencing or PCR error
  - Base call errors during sequencing could result in sequence mismatches during alignment. As these are likely to be random errors, they are expected to only be in low numbers of reads. Errors resulting in sequence shift (indels) should also be uncommon in Illumina sequencing
  - Errors during amplification steps of library preparation could result in “jackpot” mutations that propagate with further amplification cycles
2. Multiple potential alignment sites
  - The portion of a TIS read corresponding to genomic sequence is generally very short (only 16 to 17 bp for INSeq reads) so there is potential for ambiguous alignment of some sequences to more than one position
3. Reference strain vs. parent strain variation
  - Read alignment is generally performed against a completed, annotated reference sequence of the microbial strain, however there is the potential for differences between the sequence of this reference strain and the sequence of the parent strain that was used in the lab to produce the transposon mutant library. These differences could include SNVs and/or indels that could result in inaccurate alignments

  INSeq_read_preprocess.pl performs the following functions (currently only for sequencing data produced using the INSeq protocol):
  
  1. Separate multiplexed libraries by "barcode" sequence
  2. Identify and trim transposon sequence from reads
  3. Align reads to reference sequence(s)
  4. Assign aligments to positions in reference sequence
  5. Read count correction and filtering
  6. Output read counts
  
  Reads with more than one base mismatch out of 16 or 17 are removed.
  Generally, only reads aligning to potential Mariner transposon insertion sites (i.e. with terminal TA motifs) are counted unless the alignment is potentially "shifted", i.e. a perfect full-length alignment without mismatches at a potential transposon insertion site except that the alignment ends on "T" instead of "TA", or there is a single base mismatch at the terminal "T" or "A".
  Generally reads 



