# IR1-transcript-elucidation-pipeline...

A pipeline for analysing IR1 copy number and splice variation of IR1-centred transcripts in Epstein-Barr virus (EBV) long-read RNA-seq datasets. The pipeline is implemented in two R scripts that are command line executable:

1. [full_length_transcript_identification.r](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/src/full_length_long_read_identification.r) isolates full-length transcripts from long-read SAM files.

2. [IR1_transcript_variation_calculation.r](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/src/IR1_read_correction_and_elucidation.r) calculates the IR1 copy number and exon composition of full-length IR1-centred transcripts.

For exemplar data preprocessing steps see INSERT HERE.


## Full length transcript identification and isolation
