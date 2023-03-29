# IR1-transcript-elucidation-pipeline...

A pipeline for analysing IR1 copy number and splice variation of IR1-centred transcripts in Epstein-Barr virus (EBV) long-read RNA-seq datasets. The pipeline is implemented in two R scripts that are command line executable:

1. [full_length_long_read_identification.r](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/src/full_length_long_read_identification.r) isolates full-length transcripts from SAM files and outputs them in a new SAM file.

2. [IR1_read_correction_and_elucidation.r](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/src/IR1_read_correction_and_elucidation.r) calculates the IR1 copy number and exon composition of each IR1-centred transcript.
