#!/usr/bin/Rscript --vanilla

###################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# PART TWO OF TWO OF THE IR1 TRANSCRIPT ELUCIDATION PIPELINE (Aaron Logsdon, 2023. Written in R version 4.2.2 (using Bioconductor version 3.16) within Visual Studio Code version 1.76 on macOS Ventura version 13.0)
# NAME: IR1_read_correction_and_elucidation.r
# PURPOSE: reveal the exonic structure, promoter usage, 3' splice site usage, and IR1 repeat count of full-length reads from part one if used for EBV transcriptomics
# Dependencies: 
    # DEPENDENCIES: 
    # 1) R must be installed locally (see http://lib.stat.cmu.edu/R/CRAN/)
    # 2) The R packages "Biostrings", "IRanges", "stringi", and "tidyverse" must be installed using the R console (see https://www.bioconductor.org/install/)
# FUNCTIONS:
    # 1) Adds the W0 exon back to the read sequences where it has been soft-clipped, using a user-defined minimum left clip length to determine which reads to correct. We used 6 as a minimum left clip length for dRNA-seq Oxford Nanopore Technologies data
    # 2) Updates the input sam files for W0 correction and outputs (remember to add the header back) the result to the current working directory (file name format: e.g. HB9_full_length_reads_W0_corrected.sam)
    # 3) Calculates the size of each exon in IR1-containing reads for each sample using the user-defined start/end windows and splice windows and outputs the result to the current working directory (file name format: e.g. HB9_IR1_read_exons.txt). We used 30 as a start/end window and 2 as a splice window for dRNA-seq Oxford Nanopore Technologies data
    # 4) Calculates the 3' end splice site and promoter usage of IR1-containing reads
    # 4) Calculates the IR1-repeat count according to the alignment tool used by the user, for IR1-containing reads for each sample
    # 5) Calculates the IR1-repeat count according to read nucleotide count, for IR1-containing reads for each sample
    # 6) Outputs the results of the previous three steps to the current working directory (file name format: e.g. HB9_IR1_repeat_counts.txt)
# LIMITATIONS:
    # 1) The script is sepcific for EBV transcriptomic usage
    # 2) The script requires the IR1_exon_coordinates.bed file (provided in test data) to work. Exon names must be left as is but coordinates can be changed to fit the user's genome references. Also, the number of W exon groups can be adapted to match the number of IR1 repeats in the user's genome reference (ensure W0, W1, W1_prime, W2, W2_delta are present for each IR1 repeat)
    # 3) The input sam files must be headerless and headers must be added back by the user after
    # 4) Adding back W0 to reads can result in subtle variability around the W1 and W1' start sites of reads. However, this is accounted for in calculations
    # 5) The script is optimised for Oxford Nanopore Technologies and may cause suboptimal results when used with other long-read sequencing technologies
    # 6) Users are able to define start/end windows and splice windows for exon and repeat calculations. However, excessviely large windows may result in suboptimal results
# COMMAND LINE FORMAT: Rscript IR1_read_correction_and_elucidation.r sam_full_length_reads_one.sam,sam_full_length_reads_two.sam,...,sam_full_length_reads_n.sam sam_name_one,sam_name_two,...,sam_name_n sam_genome_reference_one.fasta,sam_genome_reference_two.fasta,...,sam_genome_reference_n.fasta sam_IR1_exon_coordinates_one.bed,sam_IR1_exon_coordinates_two.bed,...,sam_IR1_exon_coordinates_n.bed mininum_W0_clip_length start_end_window splice_window
# COMMAND LINE EXAMPLE: Rscript ./IR1_read_correction_and_elucidation.r ./HB9_full_length_reads.sam HB9 ./HB9_genome.fasta ./HB9_IR1_exon_coordinates.bed 6 30 2
# FOLLOW UP: independent analysis of the output files
###################################################################################################################################################################################################################################################################################################################################################################################################################################################################

### Load libraries ###
library(Biostrings)
library(IRanges)
library(stringi)
library(tidyverse)

### Initiate user input ###
args <- commandArgs(trailingOnly = TRUE)
sam_files <- strsplit(args[1], ",")[[1]]
sam_names <- strsplit(args[2], ",")[[1]]
sam_references <- strsplit(args[3], ",")[[1]]
sam_reference_IR1_coordinates <- strsplit(args[4], ",")[[1]]
mininum_W0_clip_length <- as.numeric(args[5])
start_end_window <- as.numeric(args[6])
splice_window <- as.numeric(args[7])

### Add back the small 5' W0 exon to reads from which it has been removed by soft-clipping ###
combined_reads <- vector(mode = "list", length = length(sam_files))
W0_corrected_qnames <- c()
sample_Wp_reads <- data.frame()
for(i in 1:length(sam_files)){
    sample_reads <- read.table(sam_files[i], sep = "\t", header = FALSE)
    colnames(sample_reads) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual")
    reference_IR1_coordinates <- read.table(sam_reference_IR1_coordinates[i], sep = "\t", header = TRUE)
    reference_W1_coordinates <- subset(reference_IR1_coordinates, reference_IR1_coordinates$exon == "W1")
    reference_W1_prime_coordinates <- subset(reference_IR1_coordinates, reference_IR1_coordinates$exon == "W1_prime")
    sample_Wp_reads <- data.frame(matrix(ncol = ncol(sample_reads), nrow = 0))
    for(j in 1:nrow(sample_reads)){
        for(k in 1:nrow(reference_W1_coordinates)){
            if((as.numeric(sample_reads$pos[j]) %in% (as.numeric(reference_W1_coordinates$start[k]) - start_end_window):(as.numeric(reference_W1_coordinates$start[k]) + start_end_window)) | (as.numeric(sample_reads$pos[j]) %in% (as.numeric(reference_W1_prime_coordinates$start[k]) - start_end_window):(as.numeric(reference_W1_prime_coordinates$start[k]) + start_end_window))){
                sample_Wp_reads <- rbind(sample_Wp_reads, sample_reads[j, ]) # store reads that are within the user-defined window of the W1 or W1_prime exon reference start positions
            }
        }
    }
    sample_reference <- readDNAStringSet(sam_references[i])
    sample_reference <- paste(sample_reference, collapse = "")
    sample_reference <- unlist(strsplit(as.character(sample_reference), split = "")) # split the genome reference into a vector characters
    reference_W0_coordinates <- subset(reference_IR1_coordinates, reference_IR1_coordinates$exon == "W0")
    for(j in 1:nrow(sample_Wp_reads)){
        if(sample_Wp_reads$flag[j] == "0"){ # check if the read is mapped to the plus strand
            reverse_left_clip <- ""
            left_clip_length <- 0
            left_S_index <- 0
            cigar <- unlist(strsplit(sample_Wp_reads$cigar[j], "")) # split the cigar string into a vector of characters
            for(k in 1:length(cigar)){
                if(cigar[k] == "S"){ # check for presence of left soft clip (S)
                    left_S_index <- k
                    for(m in (k - 1):1){
                        reverse_left_clip <- paste(reverse_left_clip, cigar[m], sep = "") # extract left clip numbers
                    }
                    left_clip_length <- as.numeric(stri_reverse(reverse_left_clip)) # calculate left clip length
                    break
                }
                else if(cigar[k] == "M" | cigar[k] == "I" | cigar[k] == "D" | cigar[k] == "N" | cigar[k] == "H"){ # check for absence of a left soft clip
                    break
                }
            }
            if(left_clip_length >= mininum_W0_clip_length){ # check if the left clip is long enough to be a W0 exon as defined by the user
                W0_corrected_qnames <- c(W0_corrected_qnames, sample_Wp_reads$qname[j])
                nearest_W0_index <- which.min(abs(as.numeric(reference_W0_coordinates$start) - as.numeric(sample_Wp_reads$pos[j]))) # calculate the index of the reference W0 exon that is closest to the start position of the read
                nearest_W0_sequence <- sample_reference[as.numeric(reference_W0_coordinates$start[nearest_W0_index]):as.numeric(reference_W0_coordinates$end[nearest_W0_index])]
                skip_length <- as.numeric(sample_Wp_reads$pos[j]) - as.numeric(reference_W0_coordinates$end[nearest_W0_index]) - 1 # calculate the number of bases to skip between the end of the W0 exon and the start of the read
                new_cigar <- paste(as.character(length(nearest_W0_sequence)), "M", skip_length, "N", substring(sample_Wp_reads$cigar[j], (left_S_index + 1), length(cigar)), sep = "") # create a new cigar string that includes the W0 exon
                read_sequence <- unlist(strsplit(sample_Wp_reads$seq[j], ""))
                retained_read_sequence <- substring(sample_Wp_reads$seq[j], (left_clip_length + 1), length(read_sequence))
                nearest_W0_sequence <- paste(nearest_W0_sequence, collapse = "")
                new_read_sequence <- paste(nearest_W0_sequence, retained_read_sequence, sep = "") # create a new read sequence that includes the W0 exon
                sample_reads$cigar[sample_reads$qname == sample_Wp_reads$qname[j]] <- new_cigar
                sample_reads$seq[sample_reads$qname == sample_Wp_reads$qname[j]] <- new_read_sequence
                sample_reads$pos[sample_reads$qname == sample_Wp_reads$qname[j]] <- as.numeric(reference_W0_coordinates$start[nearest_W0_index])
            }
        }
    }
    combined_reads[[i]] <- sample_reads
    write.table(combined_reads[[i]], paste0(sam_names[i], "_full_length_reads_W0_corrected.sam"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) # output the headerless sam file corrected for W0 exon clipping
}

### Calculate the size of each exon in each IR1-containing read and output the results ###
combined_IR1_read_exon_lengths <- vector(mode = "list", length = length(sam_files))
for(i in 1:length(sam_files)){
    sample_reads <- combined_reads[[i]]
    reference_IR1_coordinates <- read.table(sam_reference_IR1_coordinates[i], sep = "\t", header = TRUE)
    sample_IR1_reads <- subset(sample_reads, ((as.numeric((sample_reads$pos) %in% (as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "C1"] - start_end_window):(as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "C1"] + start_end_window))) | (as.numeric((sample_reads$pos) %in% (as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "miR-BHRF1-1_processed"] - start_end_window):(as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "miR-BHRF1-1_processed"] + start_end_window))) | (as.numeric((sample_reads$pos) %in% (as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "miR-BHRF1-2_processed"] - start_end_window):(as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "miR-BHRF1-2_processed"] + start_end_window))) | (as.numeric((sample_reads$pos) %in% (as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "miR-BHRF1-3_processed"] - start_end_window):(as.numeric(reference_IR1_coordinates$start)[reference_IR1_coordinates$exon == "miR-BHRF1-3_processed"] + start_end_window))) | sample_reads$qname %in% W0_corrected_qnames) == TRUE & sample_reads$flag == 0) # subset the sam file to include only IR1-containing reads
    sample_IR1_read_exon_lengths <- data.frame(matrix(ncol = nrow(reference_IR1_coordinates) + 2, nrow = 0))
    for(j in 1:nrow(sample_IR1_reads)){
        IR1_read_exons <- vector(mode = "character", length = nrow(reference_IR1_coordinates))
        IR1_read_exons[1:length(IR1_read_exons)] <- "0"
        cigar <- unlist(strsplit(sample_IR1_reads$cigar[j], split = "")) # split the cigar string into a vector of characters
        end_of_cigar <- FALSE
        exon_start <- as.numeric(sample_IR1_reads$pos[j])
        exon_end <- exon_start - 1
        exon_length <- 0
        first_exon_found <- FALSE
        first_exon_is_W0 <- FALSE
        letter_other_than_S_found <- FALSE
        letter_other_than_H_found <- FALSE
        for(k in 1:length(cigar)){
            reverse_event <- ""
            event_length <- 0
            if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "I" | cigar[k] == "N" | cigar[k] == "S" | cigar[k] == "H"){ # if the current character is a letter aka an event
                if(k == length(cigar)){
                    end_of_cigar <- TRUE # if the current event is the last event in the cigar string, mark the end of the cigar string
                }
                for(l in (k - 1):1){
                    if(cigar[l] != "M" & cigar[l] != "I" & cigar[l] != "D" & cigar[l] != "N" & cigar[l] != "S" & cigar[l] != "H"){ 
                        reverse_event <- paste(reverse_event, cigar[l], sep = "") # extract event numbers
                    }
                    else{
                        break
                    }
                }
                event_length <- as.numeric(stri_reverse(reverse_event)) # calculate the event length
                if(cigar[k] == "M" | cigar[k] == "D"){ # if the current event is a match or deletion
                    letter_other_than_S_found <- TRUE
                    letter_other_than_H_found <- TRUE
                    exon_end <- exon_end + event_length # add the event length to the exon end
                }
                if(cigar[k] == "M" | cigar[k] == "I" | cigar[k] == "D"){ # if the current event is a match or insertion
                    letter_other_than_S_found <- TRUE
                    letter_other_than_H_found <- TRUE
                    if(cigar[k] == "M"){
                        exon_length <- exon_length + event_length # add the event length to the exon length
                    }
                    else if((cigar[k] == "I" & event_length >= 2)){ # insertions below 2 are not added (average Oxford nanopore indel is ~1.5bp)
                        exon_length <- exon_length + event_length # add the event length to the exon length
                    }
                    else if(cigar[k] == "D" & event_length <= 2){ # deletions of 2 and less are added back (average Oxford nanopore indel is ~1.5bp)
                        exon_length <- exon_length + event_length # add the event length to the exon length
                    }

                }
                if(cigar[k] == "N" & first_exon_found == FALSE){ # if the current event is a skipped region and the first exon has not been found
                    for(l in 1:nrow(reference_IR1_coordinates)){
                        if((abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) <= start_end_window) & (abs(diff(c(exon_end, as.numeric(reference_IR1_coordinates$end[l])))) <= splice_window)){ # if the first exon is found based on a user-defined start window and splice window 
                            if(reference_IR1_coordinates$exon[l] == "W1"){
                                if(abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) < abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l + 1]))))){
                                    IR1_read_exons[l] <- exon_length # add the exon length to the exon in the IR1 read exons vector (for W1)
                                    first_exon_found <- TRUE
                                    break
                                }
                                else if(abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) > abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l + 1]))))){
                                    IR1_read_exons[l + 1] <- exon_length # add the exon length to the exon in the IR1 read exons vector (for W1_prime)
                                    first_exon_found <- TRUE
                                    break
                                }
                            }
                            IR1_read_exons[l] <- exon_length # add the exon length to the exon in the IR1 read exons vector (for any exon other than W1 or W1_prime)
                            if(reference_IR1_coordinates$exon[l] == "W0"){ 
                                first_exon_is_W0 <- TRUE # if the first exon is W0, mark that the first exon is W0
                            }
                            break
                        }
                    }
                    exon_start <- exon_end + event_length + 1 # reset the exon start and end
                    exon_end <- exon_start - 1
                    exon_length <- 0
                    if(first_exon_is_W0 == TRUE){ 
                        next # if the first exon is W0, do not mark the first exon as found
                    }
                    first_exon_found <- TRUE
                    next
                }
                else if((cigar[k] == "N" & first_exon_found == TRUE)){ # if the current event is a skipped region and the first exon has been found
                    for(l in 1:nrow(reference_IR1_coordinates)){
                        if(((abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) <= splice_window) & (abs(diff(c(exon_end, as.numeric(reference_IR1_coordinates$end[l])))) <= splice_window))){ # if the exon is found based on a user-defined splice window
                            if(reference_IR1_coordinates$exon[l] == "W1"){
                                if(abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) < abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l + 1]))))){
                                    IR1_read_exons[l] <- exon_length # add the exon length to the exon in the IR1 read exons vector (for W1)
                                    break
                                }
                                else if(abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) > abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l + 1]))))){
                                    IR1_read_exons[l + 1] <- exon_length # add the exon length to the exon in the IR1 read exons vector (for W1_prime)
                                    break
                                }
                            }
                            IR1_read_exons[l] <- exon_length # add the exon length to the exon in the IR1 read exons vector (for any exon other than W1 or W1_prime)
                            break
                        }
                    }
                    exon_start <- exon_end + event_length + 1 # reset the exon start and end
                    exon_end <- exon_start - 1
                    exon_length <- 0
                }
                else if(end_of_cigar == TRUE){ # if the end of the cigar string has been reached
                    for(l in 1:nrow(reference_IR1_coordinates)){
                        if(first_exon_found == TRUE){ # if the first exon has been found
                            if(abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) <= splice_window & abs(diff(c(exon_end, as.numeric(reference_IR1_coordinates$end[l])))) <= start_end_window){
                                IR1_read_exons[l] <- exon_length # add the exon length to the exon in the IR1 read exons vector
                                break
                            }
                        }
                        else if(first_exon_found == FALSE){ # if the first exon has not been found
                            if(abs(diff(c(exon_start, as.numeric(reference_IR1_coordinates$start[l])))) <= start_end_window & abs(diff(c(exon_end, as.numeric(reference_IR1_coordinates$end[l])))) <= start_end_window){
                                IR1_read_exons[l] <- exon_length # add the exon length to the exon in the IR1 read exons vector
                                break
                            }
                        }
                    }
                }
            }
        }
        sample_IR1_read_exon_lengths <- rbind(sample_IR1_read_exon_lengths, c(sample_IR1_reads$qname[j], IR1_read_exons))
    }
    colnames(sample_IR1_read_exon_lengths) <- c("qname", reference_IR1_coordinates$exon)
    combined_IR1_read_exon_lengths[[i]] <- sample_IR1_read_exon_lengths
    write.table(combined_IR1_read_exon_lengths[[i]], paste0(sam_names[i], "_IR1_read_exons.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

### Determine the 3' element of the IR1-containing reads ###
combined_three_prime_elements <- vector(mode = "list", length = length(sam_files))
for(i in 1:length(combined_IR1_read_exon_lengths)){
    sample_exon_lengths <- combined_IR1_read_exon_lengths[[i]]
    sample_three_prime_elements <- c()
    for(j in 1:nrow(sample_exon_lengths)){
        three_prime_element_index <- 0
        for(k in 2:ncol(sample_exon_lengths)){
            if(sample_exon_lengths[j, k] != "0"){
                three_prime_element_index <- k # if the current exon is not 0, mark the index of the current exon as the 3' element
            }
        }
        if(three_prime_element_index != 0){
            sample_three_prime_elements  <- c(sample_three_prime_elements , colnames(sample_exon_lengths)[three_prime_element_index]) # if the 3' element index is not 0, add the 3' element to the vector of 3' elements
        }
        else{
            sample_three_prime_elements  <- c(sample_three_prime_elements , NA) # if the 3' element index is 0, add NA to the vector of 3' elements
        }
    }
    combined_three_prime_elements[[i]] <- sample_three_prime_elements
}

### Calculate the number of IR1 repeats according to the alignment tool used, in IR1-containing reads ###
combined_aligner_IR1_repeat_counts <- vector(mode = "list", length = length(sam_files))
for(i in 1:length(combined_IR1_read_exon_lengths)){
    sample_W_exon_lengths <- combined_IR1_read_exon_lengths[[i]][, grepl("^W", names(combined_IR1_read_exon_lengths[[i]]))]
    sample_W_exon_lengths <- cbind(sample_W_exon_lengths, combined_IR1_read_exon_lengths[[i]]$qname)
    colnames(sample_W_exon_lengths)[ncol(sample_W_exon_lengths)] <- "qname"
    reference_IR1_coordinates <- read.table(sam_reference_IR1_coordinates[i], sep = "\t", header = TRUE)
    sample_aligner_IR1_repeat_counts <- data.frame(matrix(ncol = 4, nrow = 0))
    for(j in 1:nrow(sample_W_exon_lengths)){
        promoter <- ""
        repeat_count <- 0
        W0_count <- 0
        W0_used <- FALSE
        W1_or_W1_prime_used <- "NA"
        for(k in seq(1, ncol(sample_W_exon_lengths) - 1, by = 5)){
            
            if(sample_W_exon_lengths[j, k] == "0" & W0_used == FALSE){
                W0_count <- W0_count + 1 # increase the W0 count if the current W0 is not used in the read and W0 usage is not yet confirmed
            }
            else if(sample_W_exon_lengths[j, k] != "0" & W0_used == FALSE){
                W0_count <- W0_count + 1
                W0_used <- TRUE # confirm W0 usage if the current W0 is used in the read
            }
            if(sample_W_exon_lengths[j, k + 1] != "0" | sample_W_exon_lengths[j, k + 2] != "0"){
                repeat_count <- repeat_count + 0.5 # add 0.5 to the repeat count if the current W1 or W1' is used in the read
                if(W1_or_W1_prime_used  == "NA" & sample_W_exon_lengths[j, k + 1] != "0"){
                    W1_or_W1_prime_used <- "W1"
                }
                else if(W1_or_W1_prime_used  == "NA" & sample_W_exon_lengths[j, k + 2] != "0"){
                    W1_or_W1_prime_used <- "W1_prime"
                }
            }
            if(sample_W_exon_lengths[j, k + 3] != "0" | sample_W_exon_lengths[j, k + 4] != "0"){
                repeat_count <- repeat_count + 0.5 # add 0.5 to the repeat count if the current W2 or W2 delta is used in the read
            }
        }
        if(repeat_count == 0 & as.numeric(combined_reads[[i]]$pos[combined_reads[[i]]$qname == sample_W_exon_lengths$qname[j]]) >= as.numeric(reference_IR1_coordinates$end[reference_IR1_coordinates$exon == "Y1"])){
            repeat_count <- NA
        }
        if(W0_used == FALSE & as.numeric(combined_reads[[i]]$pos[combined_reads[[i]]$qname == sample_W_exon_lengths$qname[j]]) <= as.numeric(reference_IR1_coordinates$end[reference_IR1_coordinates$exon == "C1"])){
            promoter <- "Cp" # assign the promoter as Cp if W0 is not used and the read begins before the reference C1 exon end
        }
        else if(W0_used == FALSE & as.numeric(combined_reads[[i]]$pos[combined_reads[[i]]$qname == sample_W_exon_lengths$qname[j]]) > as.numeric(reference_IR1_coordinates$end[reference_IR1_coordinates$exon == "C1"])){
            promoter <- NA # assign the promoter as NA if W0 is not used and the read begins after the reference C1 exon end
        }
        else if(W0_used == TRUE & as.numeric(combined_reads[[i]]$pos[combined_reads[[i]]$qname == sample_W_exon_lengths$qname[j]]) > as.numeric(reference_IR1_coordinates$end[reference_IR1_coordinates$exon == "C2"])){
            promoter <- paste("Wp", W0_count, sep = "") # assign the promoter as WpX if W0 is used and the read begins after the reference C2 exon end
        }
        sample_aligner_IR1_repeat_counts <- rbind(sample_aligner_IR1_repeat_counts, c(sample_W_exon_lengths$qname[j], promoter, W1_or_W1_prime_used, repeat_count))
    }
    sample_aligner_IR1_repeat_counts <- add_column(sample_aligner_IR1_repeat_counts, combined_three_prime_elements[[i]], .after = 3)
    colnames(sample_aligner_IR1_repeat_counts) <- c("qname", "promoter", "3_prime_element", "W1_vs_W1_prime", "aligner_IR1_repeat_count")
    combined_aligner_IR1_repeat_counts[[i]] <- sample_aligner_IR1_repeat_counts
}

### Calculate the number of IR1 repeats according to nucleotide count, in IR1-containing reads ###
combined_IR1_repeat_counts <- vector(mode = "list", length = length(sam_files))
for(i in 1:length(combined_IR1_read_exon_lengths)){
    sample_IR1_read_exon_lengths <- combined_IR1_read_exon_lengths[[i]]
    reference_IR1_coordinates <- read.table(sam_reference_IR1_coordinates[[i]], sep = "\t", header = TRUE)
    sample_nucleotide_IR1_repeat_counts <- c()
    average_W1W2_length <- mean(c(198, 185.5, 193, 180.5)) # W1/W2, W1/W2_delta, W1_prime/W2, W1_prime/W2_delta
    C2_end_position <- reference_IR1_coordinates$end[reference_IR1_coordinates$exon == "C2"]
    Y1_start_position <- reference_IR1_coordinates$start[reference_IR1_coordinates$exon == "Y1"]
    for(j in 1:nrow(sample_IR1_read_exon_lengths)){
        reference_position <- as.numeric(combined_reads[[i]]$pos[combined_reads[[i]]$qname == sample_IR1_read_exon_lengths$qname[j]]) - 1
        repeat_count <- 0
        nucleotide_count <- 0
        cigar <- unlist(strsplit(combined_reads[[i]]$cigar[combined_reads[[i]]$qname == sample_IR1_read_exon_lengths$qname[j]], split = "")) # split the CIGAR string into a vector of characters
        for(k in 1:length(cigar)){
            reverse_event <- ""
            event_length <- 0
            if(reference_position <= as.numeric(C2_end_position)){ # check if the current position (on the reference genome) is not past the C2 exon end
                if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "N"){ # check for a match (M), deletion (D), or skip (N) event
                    for(l in (k - 1):1){
                        if(cigar[l] != "M" & cigar[l] != "I" & cigar[l] != "D" & cigar[l] != "N" & cigar[l] != "S" & cigar[l] != "H"){
                            reverse_event <- paste(reverse_event, cigar[l], sep = "") # extract the event numbers
                        }
                        else{
                            break
                        }
                    }
                    event_length <- as.numeric(stri_reverse(reverse_event)) # calculate the event length
                    reference_position <- as.numeric(reference_position) + event_length # update the reference position
                }
            }
            else if(reference_position > as.numeric(C2_end_position) & reference_position < (as.numeric(Y1_start_position) - 1)){ # check if the current position (on the reference genome) is past the C2 exon end
                if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "I" | cigar[k] == "N"){ # check for a match (M), deletion (D), insertion (I), or skip (N) event
                    for(l in (k - 1):1){
                        if(cigar[l] != "M" & cigar[l] != "I" & cigar[l] != "D" & cigar[l] != "N" & cigar[l] != "S" & cigar[l] != "H"){
                            reverse_event <- paste(reverse_event, cigar[l], sep = "") # extract the event numbers
                        }
                        else{
                            break
                        }
                    }
                    event_length <- as.numeric(stri_reverse(reverse_event)) # calculate the event length
                    if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "N"){ # check for a match (M), deletion (D), or skip (N) event
                        reference_position <- as.numeric(reference_position) + event_length # update the reference position
                    }
                    if(cigar[k] == "M"){ # check for a match (M) or deletion (D) event
                        if(reference_position >= (as.numeric(Y1_start_position) - 1)){ # check if the current position (on the reference genome) is past the Y1 exon start
                            nucleotide_count <- nucleotide_count + (reference_position - (as.numeric(Y1_start_position) - 1)) # add the correct number of nucleotides to the count
                            break
                        }
                        nucleotide_count <- nucleotide_count + event_length # add the event length to the count
                    }
                    else if(cigar[k] == "I"){ # insertions below 2 are not added (average Oxford nanopore indel is ~1.5bp)
                        if(event_length >= 2){ # check if the insertion is of significant length to add (setting as 10 removes sequencing error noise and allows us to concentrate on interesting events like extra W exons)
                            if(reference_position >= (as.numeric(Y1_start_position) - 1)){ # check if the current position (on the reference genome) is past the Y1 exon start
                                nucleotide_count <- nucleotide_count + (reference_position - (as.numeric(Y1_start_position) - 1)) # add the correct number of nucleotides to the count
                                break
                            }
                            nucleotide_count <- nucleotide_count + event_length # add the event length to the count
                        }
                    }
                    else if(cigar[k] == "D"){
                        if(event_length <= 2){ # deletions of 2 and below are added back (average Oxford nanopore indel is ~1.5bp)
                            if(reference_position >= (as.numeric(Y1_start_position) - 1)){ # check if the current position (on the reference genome) is past the Y1 exon start
                                nucleotide_count <- nucleotide_count + (reference_position - (as.numeric(Y1_start_position) - 1)) # add the correct number of nucleotides to the count
                                break
                            }
                            nucleotide_count <- nucleotide_count + event_length
                        }
                    }
                    else if(cigar[k] == "N"){ # check for a skip (N) event
                        if(reference_position >= (as.numeric(Y1_start_position) - 1)){ # check if the current position (on the reference genome) is past the Y1 exon start
                            break
                        }
                    }
                }
            }
        }
        if(sample_IR1_read_exon_lengths$qname[j] %in% W0_corrected_qnames){
            nucleotide_count <- nucleotide_count - 27 # adjust nucleotide count for W0 in relevant reads
        }
        if(nucleotide_count > 0){ # check if the nucleotide count is greater than 0
            repeat_count <- round(as.numeric(nucleotide_count) / average_W1W2_length, digits = 1) # calculate the repeat count
        }
        else if(nucleotide_count == 0){ # check if the nucleotide count is equal to 0
            repeat_count <- NA
        }
        sample_nucleotide_IR1_repeat_counts <- c(sample_nucleotide_IR1_repeat_counts, repeat_count)
    }
    combined_IR1_repeat_counts[[i]] <- cbind(combined_aligner_IR1_repeat_counts[[i]], sample_nucleotide_IR1_repeat_counts)
    colnames(combined_IR1_repeat_counts[[i]]) <- c(colnames(combined_aligner_IR1_repeat_counts[[i]]), "nucleotide_IR1_repeat_count")
    write.table(combined_IR1_repeat_counts[[i]], file = paste0(sam_names[i], "_IR1_repeat_counts.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}