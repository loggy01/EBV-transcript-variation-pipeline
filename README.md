# EBV transcript variation pipeline
A pipeline for analysing the diversity of IR1 copy number and alternative splicing of Cp- and Wp-initiated transcripts in Epstein-Barr virus (EBV) long-read RNA-seq datasets. The pipeline is implemented in two R scripts which can be downloaded and ran on the command line:
1. [full_length_transcript_identification.r](https://github.com/loggy01/IR1-transcript-variation-pipeline/blob/main/src/full_length_transcript_identification.r) isolates and extracts full-length transcripts from EBV long-read RNA-seq SAM files.
2. [transcript_variation_calculation.r](https://github.com/loggy01/IR1-transcript-variation-pipeline/blob/main/src/transcript_variation_calculation.r) calculates the diversity of IR1 copy number and exon content of each extracted full-length Cp- and Wp-initiated transcript.

For necessary data preprocessing steps see [here](https://github.com/loggy01/IR1-transcript-variation-pipeline/blob/main/Additional%20files/command_lines.docx) and [here](https://github.com/loggy01/IR1-transcript-variation-pipeline/blob/main/Additional%20files/bam_filtration.r).


## full_length_transcript_identification.r

### Purpose
To separate each of one or more EBV long-read RNA-seq headerless SAM files into four output SAM files with the following content respectively:
1. Full length transcripts
2. 3' intact transcripts
3. 5' intact transcripts
4. Unassigned reads

### Dependencies
[R](http://lib.stat.cmu.edu/R/CRAN/) must be installed locally along with the packages [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [stringi](https://cran.r-project.org/web/packages/stringi/index.html), and [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html).

### Functions
1. Capture stage one: defines each 5' and 3' end as reproducible if they occur in a (user-defined) large enough cluster in a (user-defined) small enough window size.
2. Capture stage two: redefines read ends labelled as noise in the previous stage as reproducible if they occur within a (user-defined) small enough window size of a reproducible cluster.  
3. Capture stage three: adds back soft-clipped bases to read ends still defined as noise and redfines them as reproducible if they now occur within a (user-defined) small enough window size of a reproducible cluster. 
4. Assigns reads to one of four output SAM files based on whether a read has a reproducile 5' and/or 3' end.

### Limitations
1. Secondary and supplementary reads are not supported and must be removed from SAM files prior. See [here](https://github.com/loggy01/IR1-transcript-variation-pipeline/blob/main/Additional%20files/command_lines.docx).
2. Padded sequences are not supported.
3. SAM files must have cigar strings in original format (M, not = and X).
4. Each input SAM file must have its header removed prior. See [here](https://github.com/loggy01/IR1-transcript-variation-pipeline/blob/main/Additional%20files/command_lines.docx).

### Input and output
Input: 
1. One or more SAM files (list of directories).
2. Matching list of SAM file names (list of strings). 
3. Window size for capture stage one (integer). We use 2.
4. Minimum cluster size for capture stage one (integer). We use 5.
6. Window size for capture stage two (integer). We use 20.
7. Window size for capture stage three (integer). We use 20.

Output (per input headerless SAM file):
1. ./sam_name_full_length_transcripts.sam
2. ./sam_name_3'_intact_transcripts.sam
3. ./sam_name_5'_intact_transcripts.sam
4. ./sam_name_unassigned_reads.sam

### Command line
````shell
Rscript ./full_length_transcript_identification.r ./sam_one.sam,./sam_two.sam name_one,name_two 2 5 20 20
````


## transcript_variation_calculation.r

### Purpose
To calculate IR1 copy number and alternative splicing diversity of Cp- and Wp-initiated transcripts in full-length transcript headerless SAM files

### Dependencies
[R](http://lib.stat.cmu.edu/R/CRAN/) must be installed locally along with [Bioconductor](https://www.bioconductor.org/install/) and the packages [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html), [plyr](https://cran.r-project.org/web/packages/plyr/index.html), [stringi](https://cran.r-project.org/web/packages/stringi/index.html), and [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html).

### Functions
1. Adds back the W0 exon to transcripts in which it has been erroneously clipped off, using a (user-defined) mininum soft-clip size for consideration, and outputs the updated SAM files.
2. Calculates the aligned exon content in each Cp- and Wp-initiated transcript using two (user-defined) window sizes to account for position variability of 5' and 3' exons and internal exons respectively. 
3. Calculates the diversity of IR1 copy number in Cp- and Wp-initiated transcripts based on the number of aligned W1, W1', W2, and W2' exons for each read. This is termed the *alignment count* and only includes aligned IR1 copies.
4. Calculates the diversity of IR1 copy number in Cp- and Wp-initiated transcripts based on the size of the aligned IR1 compoent of each read. This is termed the *distance count* and includes aligned and misaligned IR1 copies.
5. Outputs a text file containing the aligned exon content of each read and another text file containing IR1 copy number diversity as calculated by the alignment an distance count, for each input headerless SAM file.

### Limitations 
1. The script requires the pipeline's EBV exon coordinates BED file to work (see [here](https://github.com/loggy01/EBV-transcript-variation-pipeline/blob/main/examples/transcript_variation_calculation/input.bed)), which cannot be changed a part from to add or remove IR1 copies (each copy must have an entry for W0, W1, W1_prime, W2, W2_prime). The default is six copies.
2. The script is optimised for Oxford Nanopore Technologies data and may perform sub-optimally with other long-read RNA-seq technologies.

### Input and output
Input: 
1. One or more full-length transcript SAM files from full_length_transcript_identification.r (list of directories)
2. Matching list of SAM file name(s) (list of strings).
3. Matching list of reference genome FASTA file(s) (list of directories)
4. Matching list of EBV exon coordinates BED file(s) (list of directories)
5. Minimum clip size for W0 add back consideration (integer). We use 5.
6. Window size for 5' and 3' exon position variability from the associated EBV exon cooridnates BED file (integer). We use 20.
7. Window size for internal exon position variability from the associated EBV exon cooridnates BED file (integer). We use 2.

Output (per input headerless SAM file):
1. ./sample_name_full_length_transcripts_final.sam (remember to add back the header yourself)
2. ./sam_name_IR1_exon_compositions.txt
3. ./sam_name_IR1_copy_number_variation.txt

### Command line
````shell
Rscript ./transcript_variation_calculation.r ./full_length_transcripts_one.sam,./full_length_transcripts_two.sam name_one,name_two ./ref_one.fa,./ref_two.fa ./coordinates_one.bed,./coordinates_two.bed 5 20 2
