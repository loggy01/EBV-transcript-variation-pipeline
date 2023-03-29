# IR1 transcript variation pipeline
A pipeline for analysing IR1 copy number and splice variation of IR1-centred transcripts in Epstein-Barr virus long-read RNA-seq datasets. The pipeline is implemented in two R scripts which can be downloaded and run on the command line:
1. [full_length_transcript_identification.r](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/src/full_length_long_read_identification.r) isolates full-length transcripts from long-read SAM files.
2. [IR1_transcript_variation_calculation.r](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/src/IR1_read_correction_and_elucidation.r) calculates the IR1 copy number and exon composition of full-length IR1-centred transcripts.

For exemplar data preprocessing steps see [here]() and [here]().


## full_length_transcript_identification.r

### Purpose
To separate each of an *n* number of long-read RNA-seq SAM files into four output SAM files:
1. Full length transcripts
2. 3' intact transcripts
3. 5' intact transcripts
4. Unassigned reads

### Dependencies
[R](http://lib.stat.cmu.edu/R/CRAN/) must be installed locally along with the packages [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [stringi](https://cran.r-project.org/web/packages/stringi/index.html), and [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html).

### Functions
1. Capture stage one: captures reproducible 5' and 3' ends based on a user-defined mininum group size occuring within a user-defined window size (±).
3. Capture stage two: captures uncaptured 5' and 3' ends within a user-defined window size (±) from the centre of a reproducible group.
4. Capture stage three: reverses soft-clipping and captures 5' and 3' ends now within a user-defined window size (±).
5. Assigns reads to one of four output SAM files based on whether a read has a reproducile 5' and/or 3' end.

### Limitations
1. Secondary and supplementary reads are not supported and must be removed from SAM files prior. See [here]().
2. Padded sequences are not supported.
3. SAM files must have cigar strings in original format (M, not = and X).
4. Each input SAM file must have its header removed prior. See [here]().

### Input and output
Input: 
1. *n* number of SAM files (list of directories).
2. Matching list of SAM file names (list of strings). 
3. Minimum group size for capture stage one (integer). We use 5.
4. Window size for capture stage one (integer). We use 2.
5. Window size for capture two (integer). We use 20.
6. Window size for capture stage three (integer). We use 20.

Output (four text files (for each input SAM) in your working directory):
1. ./sam_name_full_length_transcripts.sam
2. ./sam_name_3_prime_intact_transcripts.sam
3. ./sam_name_5_prime_intact_transcripts.sam
4. ./sam_name_unassigned_reads.sam

### Command line
````shell
Rscript ./full_length_transcript_identification.r ./sam_one.sam,./sam_two.sam name_one,name_two 2 5 20 20
````


## IR1_transcript_variation_calculation.r

### Purpose
To calculate IR1 copy number and splice variation of full-length IR1-centred transcripts from ./full_length_transcripts.sam.

### Dependencies
[R](http://lib.stat.cmu.edu/R/CRAN/) must be installed locally along with [Bioconductor](https://www.bioconductor.org/install/) and the packages [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html), [stringi](https://cran.r-project.org/web/packages/stringi/index.html), and [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)

### Functions
1. Adds back the W0 exon to transcripts in which it has been erroneously clipped off, using a user-defined mininum soft-clip size for consideration, and      outputs the updated SAM files.
2. Calculates absence/presence of each aligned exon in each IR1-centred transcript, using user-defined start, end, and splice windows (±).
3. Calculates IR1 copy number based on a transcript's aligned exons (alignment count).
4. Calculates IR1 copy number based on distance between the start and the end of the IR1 component of a transcript (distance count).
5. Outputs IR1 copy number variation text file and exon composition text file for each input SAM file.

### Limitations 
1. The script requires the pipeline's IR1 exon coordinates BED file to work (see [here](https://github.com/loggy01/IR1-transcript-elucidation-pipeline/blob/main/examples/IR1_read_correction_and_elucidation/input.bed), which cannot be changed a part from to add or remove IR1 copies (each copy must have an entry for W0, W1, W1_prime, W2, W2_prime).
2. The script is optimised for Oxford Nanopore Technologies data and may perform sub-optimally with other long-read RNA-seq technologies.

### Input and output
Input: 
1. *n* number of ./full_length_transcripts.sam from full_length_transcript_identification.r (list of directories)
2. Matching list of SAM file names (list of strings).
3. Matching list of reference genome FASTA files (list of directories)
4. Matching list of IR1 exon coordinates BED file (list of directories)
5. Minimum clip size for W0 add back consideration (integer). We use 5.
6. Window size for transcript start and end variability from matching IR1 cooridnates exon BED file (integer). We use 20.
7. Window size for exon splice position variability from matching IR1 cooridnates exon BED file (integer). We use 2.

Output (three text files (for each input SAM) in your working directory):
1. ./sample_name_full_length_transcripts_final.sam (remember to add back the header yourself)
2. ./sam_name_IR1_exon_compositions.txt
3. ./sam_name_IR1_copy_number_variation.txt

### Command line
````shell
Rscript ./IR1_transcript_variation_calculation.r ./full_length_transcripts_one.sam,./full_length_transcripts_two.sam name_one,name_two ./ref_one.fa,./ref_two.fa ./coordinates_one.bed,./coordinates_two.bed 5 20 2
