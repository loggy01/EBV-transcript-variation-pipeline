#!/usr/bin/Rscript --vanilla

#################################################################################################################################################################################################################################################
# PART ONE OF TWO OF THE IR1 TRANSCRIPT ELUCIDATION PIPELINE (Aaron Logsdon, 2023. Written in R version 4.2.2 using Bioconductor version 3.16 within Visual Studio Code version 1.76 on macOS Ventura version 13.0)
# NAME: full_length_long_read_identification.r
# PURPOSE: separate n number of sam files containing long read data into four read sets for each sam file:
    # 1) full length reads
    # 2) 5' end broken reads
    # 3) 3' end broken reads
    # 4) unassigned reads
# DEPENDENCIES: 
    # 1) R must be installed locally (see http://lib.stat.cmu.edu/R/CRAN/)
    # 2) The R packages "data.table", "stringi", and "tidyverse" must be installed using the R console (see https://www.bioconductor.org/install/)
# FUNCTIONS: 
    # 1) Extracts read start and end positions from n number of input sam files
    # 2) Assigns significant windows of read starts/ends using a user-defined window size (± from reference position) and user-defined minimum window read count. We used a window size of 2 and minimum window read count of 5
    # 3) Assigns unassigned and potentially significant read starts/ends missed by the previous step to significant windows based on a user defined window size (± from center of significant window in question). We used a window size of 10
    # 4) Assigns read starts/ends missed by the previous steps due to soft clipping to significant windows, again using a user defined window size (± from center of significant window in question). We used a window size of 10
    # 5) Updates each input sam file for soft clipping correction from the previous step
    # 6) Assigns reads to the mentioned read sets based on whether they appear in significant start and end windows, only significant start windows, only significant end windows, or neither respectively
    # 7) Outputs the four read sets from each sam file into the current working directory (file name format: e.g. HB9_full_length_reads.sam)
# LIMITATIONS:
    # 1) Secondary and supplementary alignments are not supported and must be removed prior using samtools on the command line (samtools view -F 256 -F 2048 input_sam > output_sam)
    # 2) Padding of sequences is not supported (this is only relevant if multiple sequence alignment was performed prior)
    # 3) Cigar strings must be in the original format aka using M instead of = and X (most alignment tools allow the choice of format)
    # 4) Each input sam file must be edited prior to remove the qual column (the script adds this back) and header (the user must add this back after) using samtools on the command line (samtools view input_sam | cut -f 1-10 > output_sam)
# COMMAND LINE FORMAT: Rscript full_length_long_read_identification.r sam_one.sam,sam_two.sam,...,sam_n.sam sam_name_one,sam_name_two,...,sam_name_n significance_window mininum_significant_read_count noise_window clip_window
# COMMAND LINE EXAMPLE: Rscript ./full_length_long_read_identification.r ./HB9.sam HB9 2 5 10 10
# FOLLOW UP: Output full length read sam files should be passed to the next script in the pipeline (IR1_transcript_correction_and_elucidation.r)
#################################################################################################################################################################################################################################################

### Load packages ###
library(data.table)
library(stringi)
library(tidyverse)

### Initiate user input ###
args <- commandArgs(trailingOnly = TRUE)
sam_files <- strsplit(args[1], ",")[[1]]
sam_names <- strsplit(args[2], ",")[[1]]
significance_window <- as.numeric(args[3])
mininum_significant_read_count <- as.numeric(args[4])
noise_window <- as.numeric(args[5])
clip_window <- as.numeric(args[6])

### Identify and store read starts (plus strand) and ends (minus strand) ###
position_data <- data.frame(matrix(nrow = 0, ncol = 6))
for(i in 1:length(sam_files)){
    sam_file <- read.table(sam_files[i], sep = "\t", header = FALSE)
    colnames(sam_file) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq")
    sam_name <- sam_names[i]
    for(j in 1:nrow(sam_file)){
        position_data <- rbind(position_data, c(sam_file$qname[j], sam_name, sam_file$flag[j], sam_file$pos[j], sam_file$cigar[j], as.character(sam_file$seq[j]))) # extracting read relevant information from sam file
    }
}
colnames(position_data) <- c("qname", "sample", "dna_strand", "bam_read_position", "read_cigar", "read_sequence")
position_data$dna_strand[position_data$dna_strand == 0] <- "+"
position_data$dna_strand[position_data$dna_strand == 16] <- "-"

### Identify read starts (minus strand) and ends (plus strand) ###
reverse_bam_read_position <- c()
for(i in 1:nrow(position_data)){
    length_on_genome <- 0
    cigar <- unlist(strsplit(position_data$read_cigar[i], split = "")) # split cigar string into characters
    for(j in 1:length(cigar)){
        if(cigar[j] == "N" | cigar[j] == "M" | cigar[j] == "D"){ # Look for skip (N), match (M), or deletion (D) events
            reverse_event_length <- ""
            event_length <- 0
            for(k in (j-1):1){
                if(cigar[k] != "S" & cigar[k] != "H" & cigar[k] != "M" & cigar[k] != "I" & cigar[k] != "D" & cigar[k] != "N"){
                    reverse_event_length <- paste(reverse_event_length, cigar[k], sep = "") # extract event numbers
                }
                else{
                    break
                }
            }
            event_length <- as.numeric(stri_reverse(reverse_event_length)) # calculate event length from numbers
            length_on_genome <- length_on_genome + event_length # track current read position on genome
        }
    }
    new_read_position <- (length_on_genome - 1) + as.numeric(position_data$bam_read_position[i])
    reverse_bam_read_position <- c(reverse_bam_read_position, new_read_position) # store the newly calculated read position on genome
}
position_data <- add_column(position_data, reverse_bam_read_position, .after = 4)

### Separate + strand and - strand reads then separate into read starts and ends ###
plus_position_data <- subset(position_data, position_data$dna_strand == "+", select = -dna_strand)
colnames(plus_position_data) <- c("qname", "sample", "read_start", "read_end", "read_cigar", "read_sequence")
minus_position_data <- subset(position_data, position_data$dna_strand == "-", select = -dna_strand)
colnames(minus_position_data) <- c("qname", "sample", "read_end", "read_start", "read_cigar", "read_sequence")
minus_position_data <- minus_position_data[, colnames(plus_position_data)] # reorder columns to match plus_position_data
combined_position_data <- list(plus_start_position_data = subset(plus_position_data, select = -read_end), plus_end_position_data = subset(plus_position_data, select = -read_start), minus_start_position_data = subset(minus_position_data, select = -read_end), minus_end_position_data = subset(minus_position_data, select = -read_start))
combined_position_data[[1]] <- combined_position_data[[1]][order(as.numeric(combined_position_data[[1]]$read_start)), ]
combined_position_data[[2]] <- combined_position_data[[2]][order(as.numeric(combined_position_data[[2]]$read_end)), ]
combined_position_data[[3]] <- combined_position_data[[3]][order(as.numeric(combined_position_data[[3]]$read_start)), ]
combined_position_data[[4]] <- combined_position_data[[4]][order(as.numeric(combined_position_data[[4]]$read_end)), ]

### Assign significant windows of read starts/ends using a user-defined significance window size and count ###
combined_position_comparison_data <- vector(mode = "list", length = 4)
unassigned_qnames <- list(plus_start = position_data$qname, plus_end = position_data$qname, minus_start = position_data$qname, minus_end = position_data$qname)
for(i in 1:length(combined_position_data)){
    comparison_data <- data.frame(matrix(nrow = 0, ncol = 3))
    last_query_index <- 0
    for(j in 1:nrow(combined_position_data[[i]])){
        significant_read_count <- 0
        significant_read_qnames <- c()
        if(j == last_query_index + 1){
            reference_position <- as.numeric(combined_position_data[[i]][j, 3]) # define reference position
            for(k in 1:nrow(combined_position_data[[i]])){
                query_position <- as.numeric(combined_position_data[[i]][k, 3]) # define query position
                if((query_position >= reference_position - significance_window) & (query_position <= reference_position + significance_window)){ # check if query position is equal to reference position ± significance_window
                    significant_read_count <- significant_read_count + 1
                    significant_read_qnames <- paste(significant_read_qnames, combined_position_data[[i]]$qname[k], sep = ",")
                    last_query_index <- k # keep track of last query index to avoid redundant comparisons
                }
            }
            if(significant_read_count >= mininum_significant_read_count){ # check if read count is greater than or equal to user_defined mininum ciunt before adding to comparison_data
                comparison_data <- rbind(comparison_data, c(reference_position, significant_read_count, significant_read_qnames))
                significant_read_qnames <- unlist(strsplit(significant_read_qnames, split = ","))
                unassigned_qnames[[i]] <- unassigned_qnames[[i]][!(unassigned_qnames[[i]] %in% significant_read_qnames)] # remove the newly assigned significant qnames from the unassigned_qnames for the current strand (plus/minus) and read type (start/end)
            }
        }
    }
    comparison_data <- cbind(comparison_data, noisy_read_count = 0, noisy_read_qnames = "", clipped_read_count = 0, clipped_read_qnames = "") # add significant data just calculated and empty columns for noisy and clipped read correction later
    colnames(comparison_data) <- c("genomic_position", "significant_read_count", "significant_read_qnames", "noisy_read_count", "noisy_read_qnames", "clipped_read_count", "clipped_read_qnames")
    comparison_data <- comparison_data[order(as.numeric(comparison_data$genomic_position)), ]
    combined_position_comparison_data[[i]] <- comparison_data
}

### Capture significant "noise" around the assigned significant windows using a user-defined window size ###
for(i in 1:length(combined_position_data)){
    assigned_noisy_qnames <- c()
    for(j in 1:nrow(combined_position_data[[i]])){
        for(k in 1:length(unassigned_qnames[[i]])){
            if(unassigned_qnames[[i]][k] == combined_position_data[[i]]$qname[j]){ 
                for(l in 1:nrow(combined_position_comparison_data[[i]])){
                    if((as.numeric(combined_position_data[[i]][j, 3]) >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[l]) - noise_window) & (as.numeric(combined_position_data[[i]][j, 3]) <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[l]) + noise_window)){ # check if noisy read is equal to an assigned significant read position ± noise_window
                        combined_position_comparison_data[[i]]$noisy_read_count[l] <- as.numeric(combined_position_comparison_data[[i]]$noisy_read_count[l]) + 1
                        combined_position_comparison_data[[i]]$noisy_read_qnames[l] <- paste(combined_position_comparison_data[[i]]$noisy_read_qnames[l], unassigned_qnames[[i]][k], sep = ",")
                        assigned_noisy_qnames <- c(assigned_noisy_qnames, unassigned_qnames[[i]][k])
                    }
                }
            }
        }
    }
    if(length(assigned_noisy_qnames) > 0){
        unassigned_qnames[[i]] <- unassigned_qnames[[i]][!(unassigned_qnames[[i]] %in% assigned_noisy_qnames)] # remove the newly assigned noisy qnames from the unassigned_qnames for the current strand (plus/minus) and read type (start/end)
    }
    combined_position_comparison_data[[i]]$noisy_read_qnames[combined_position_comparison_data[[i]]$noisy_read_qnames == ""] <- "NA"
    combined_position_comparison_data[[i]] <- combined_position_comparison_data[[i]][order(as.numeric(combined_position_comparison_data[[i]]$genomic_position)), ]
}

### Capture significant clipped reads with clips falling around significant windows, using a user-defined window size, and write read counts from TASK FOUR, TASK FIVE, and TASK SIX. Also update the input sam files for the clip corrections ###
combined_sam_files <- vector(mode = "list", length = length(sam_files))
for(i in 1:length(sam_files)){
    combined_sam_files[[i]] <- read.table(sam_files[i], sep = "\t", header = FALSE) # store the input sam files in a list
    colnames(combined_sam_files[[i]]) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq")
}
for(i in 1:length(combined_position_data)){
    assigned_clipped_qnames <- c()
    for(j in 1:nrow(combined_position_data[[i]])){
        for(k in 1:length(unassigned_qnames[[i]])){
            if(unassigned_qnames[[i]][k] == combined_position_data[[i]]$qname[j]){
                reverse_left_clip <- ""
                reverse_right_clip <- ""
                left_S_index = 0
                right_S_index = 0
                left_clip_length <- 0
                right_clip_length <- 0
                new_start_position <- 0
                new_end_position <- 0
                next_clip_is_right <- FALSE
                cigar <- unlist(strsplit(combined_position_data[[i]]$read_cigar[j], split = "")) # split cigar string into individual characters
                for(l in 1:length(cigar)){
                    if(i == 1 | i == 4){ # Look for left clipping (1 = plus strand start, 4 = minus strand end)
                        if(cigar[l] == "S"){ # check for presence of left soft clip (S)
                            left_S_index <- l 
                            for(m in (l - 1):1){
                                reverse_left_clip <- paste(reverse_left_clip, cigar[m], sep = "") # extract left clip numbers
                            }
                            left_clip_length <- as.numeric(stri_reverse(reverse_left_clip)) # calculate left clip length
                            if(i == 1){
                                new_start_position <- as.numeric(combined_position_data[[i]][j, 3]) - left_clip_length # calculate soft clip corrected start position for plus strand start
                            }
                            else if(i == 4){
                                new_end_position <- as.numeric(combined_position_data[[i]][j, 3]) - left_clip_length # calculate soft clip corrected end position for minus strand end
                            }
                            for(m in 1:nrow(combined_position_comparison_data[[i]])){
                                if(i == 1){
                                    if((new_start_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_start_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){ # check if left clip corrected plus strand start read position is equal to an assigned significant start position ± clip_window
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_start_position) # update the bam_read_position column in the position_data data frame with the soft clip corrected start position
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(as.character(left_clip_length), "M", substring(combined_position_data[[i]]$read_cigar[j], (left_S_index + 1), length(cigar)), sep = "") # update the read_cigar column in the position_data data frame with the left soft clip corrected cigar string
                                        for(n in 1:length(combined_sam_files)){
                                            if(combined_position_data[[i]]$qname[j] %in% combined_sam_files[[i]]$qname){
                                                combined_sam_files[[i]]$pos[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_start_position) # update the pos column in the sam_file data frame with the soft clip corrected start position
                                                combined_sam_files[[i]]$cigar[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- paste(as.character(left_clip_length), "M", substring(combined_position_data[[i]]$read_cigar[j], (left_S_index + 1), length(cigar)), sep = "") # update the cigar column in the sam_file data frame with the left soft clip corrected cigar string
                                            }
                                        }
                                    }
                                }
                                else if(i == 4){
                                    if((new_end_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_end_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){ # check if left clip corrected minus strand end read position is equal to an assigned significant end position ± clip_window
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_end_position) # update the bam_read_position column in the position_data data frame with the soft clip corrected end position
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(as.character(left_clip_length), "M", substring(combined_position_data[[i]]$read_cigar[j], (left_S_index + 1), length(cigar)), sep = "") # update the read_cigar column in the position_data data frame with the left soft clip corrected cigar string
                                        for(n in 1:length(combined_sam_files)){
                                            if(combined_position_data[[i]]$qname[j] %in% combined_sam_files[[i]]$qname){
                                                combined_sam_files[[i]]$pos[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_end_position) # update the pos column in the sam_file data frame with the soft clip corrected end position
                                                combined_sam_files[[i]]$cigar[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- paste(as.character(left_clip_length), "M", substring(combined_position_data[[i]]$read_cigar[j], (left_S_index + 1), length(cigar)), sep = "") # update the cigar column in the sam_file data frame with the left soft clip corrected cigar string
                                            }
                                        }
                                    }
                                }
                            }
                            break
                        }
                        else if(cigar[l] == "M" | cigar[l] == "I" | cigar[l] == "D" | cigar[l] == "N" | cigar[l] == "H"){ # check for absence of a left soft clip
                            break
                        }
                    }
                    else if(i == 2 | i == 3){ # Look for right clipping (2 = plus strand end, 3 = minus strand start)
                        if(cigar[l] == "M" | cigar[l] == "I" | cigar[l] == "D" | cigar[l] == "N" | cigar[l] == "H"){ # check for absence of a right soft clip
                            next_clip_is_right <- TRUE
                        }
                        else if(cigar[l] == "S" & next_clip_is_right == TRUE){ # check for presence of right soft clip (S)
                            right_S_index <- l
                            for(m in (l -1):1){
                                if(cigar[m] != "M" & cigar[m] != "I" & cigar[m] != "D" & cigar[m] != "N" & cigar[m] != "S" & cigar[m] != "H"){
                                    reverse_right_clip <- paste(reverse_right_clip, cigar[m], sep = "") # extract right clip numbers
                                }
                                else{
                                    break
                                }
                            }
                            right_clip_length <- as.numeric(stri_reverse(reverse_right_clip)) # calculate right clip length
                            if(i == 2){
                                new_end_position <- as.numeric(combined_position_data[[i]][j, 3]) + right_clip_length # calculate soft clip corrected end position for plus strand end
                            }
                            else if(i == 3){
                                new_start_position <- as.numeric(combined_position_data[[i]][j, 3]) + right_clip_length # calculate soft clip corrected start position for minus strand start
                            }
                            for(m in 1:nrow(combined_position_comparison_data[[i]])){
                                if(i == 2){
                                    if((new_end_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_end_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){ # check if right clip corrected plus strand end read position is equal to an assigned significant start position ± clip_window
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$reverse_bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_end_position) # update the bam_read_position column in the position_data data frame with the soft clip corrected end position
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(substring(combined_position_data[[i]]$read_cigar[j], 1, (right_S_index - 1)), "M", sep = "") # update the read_cigar column in the position_data data frame with the right soft clip corrected cigar string
                                        for(n in 1:length(combined_sam_files)){
                                            if(combined_position_data[[i]]$qname[j] %in% combined_sam_files[[i]]$qname){
                                                combined_sam_files[[i]]$pos[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_end_position) # update the pos column in the sam_file data frame with the soft clip corrected end position
                                                combined_sam_files[[i]]$cigar[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- paste(substring(combined_position_data[[i]]$read_cigar[j], 1, (right_S_index - 1)), "M", sep = "") # update the cigar column in the sam_file data frame with the right soft clip corrected cigar string
                                            }
                                        }
                                    }
                                }
                                else if(i == 3){
                                    if((new_start_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_start_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){ # check if right clip corrected minus strand start read position is equal to an assigned significant end position ± clip_window
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$reverse_bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_start_position) # update the bam_read_position column in the position_data data frame with the soft clip corrected start position
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(substring(combined_position_data[[i]]$read_cigar[j], 1, (right_S_index - 1)), "M", sep = "") # update the read_cigar column in the position_data data frame with the right soft clip corrected cigar string
                                        for(n in 1:length(combined_sam_files[[i]])){
                                            if(combined_position_data[[i]]$qname[j] %in% combined_sam_files[[i]]$qname){
                                                combined_sam_files[[i]]$pos[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_start_position) # update the pos column in the sam_file data frame with the soft clip corrected start position
                                                combined_sam_files[[i]]$cigar[combined_sam_files[[i]]$qname == combined_position_data[[i]]$qname[j]] <- paste(substring(combined_position_data[[i]]$read_cigar[j], 1, (right_S_index - 1)), "M", sep = "") # update the cigar column in the sam_file data frame with the right soft clip corrected cigar string
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(length(assigned_clipped_qnames[[i]]) > 0){
        unassigned_qnames[[i]] <- unassigned_qnames[[i]][!(unassigned_qnames[[i]] %in% assigned_clipped_qnames)] # remove newly assigned clipped reads from unassigned_qnames
    }
    combined_position_comparison_data[[i]]$clipped_read_qnames[combined_position_comparison_data[[i]]$clipped_read_qnames == ""] <- "NA"
    combined_position_comparison_data[[i]] <- combined_position_comparison_data[[i]][order(as.numeric(combined_position_comparison_data[[i]]$genomic_position)), ]
}

### Separate reads into a dual assigned significant read set, start only, end only, and an unassigned read set ###
combined_assigned_qnames <- list(plus_start_assigned_qnames = c(), plus_end_assigned_qnames = c(), minus_start_assigned_qnames = c(), minus_end_assigned_qnames = c())
plus_strand_read_sets <- vector(mode = "list", length = 4)
minus_strand_read_sets <- vector(mode = "list", length = 4)
for(i in 1:length(combined_position_comparison_data)){
    for(j in 1:nrow(combined_position_comparison_data[[i]])){
        significant_read_count <- unlist(strsplit(combined_position_comparison_data[[i]]$significant_read_qnames[j], ","))
        noisy_read_count <- unlist(strsplit(combined_position_comparison_data[[i]]$noisy_read_qnames[j], ","))
        clipped_read_count <- unlist(strsplit(combined_position_comparison_data[[i]]$clipped_read_qnames[j], ","))
        combined_assigned_qnames[[i]] <- c(combined_assigned_qnames[[i]], significant_read_count, noisy_read_count, clipped_read_count)
    }
    combined_assigned_qnames[[i]] <- combined_assigned_qnames[[i]][!combined_assigned_qnames[[i]] == "NA" & !combined_assigned_qnames[[i]] == ""]
    if(i == 2){ # Checking for completion of plus strand read assignment data collection
        plus_qname_identification_data <- rbindlist(combined_position_data, idcol = "Identification_count", fill = TRUE)[qname %in% combined_assigned_qnames[[1]] & qname %in% combined_assigned_qnames[[2]]] # isolate read qnames that have been assigned both a significant start and significant end
        plus_dual_assigned_qnames <- unique(plus_qname_identification_data$qname)
        plus_strand_read_sets[[1]] <- subset(position_data, qname %in% plus_dual_assigned_qnames) # store dual assigned reads (significant start and end)
        plus_strand_read_sets[[2]] <- subset(position_data, qname %in% combined_assigned_qnames[[1]] & !(qname %in% plus_dual_assigned_qnames)) # store start only assigned reads (insignificant end)
        plus_strand_read_sets[[3]] <- subset(position_data, qname %in% combined_assigned_qnames[[2]] & !(qname %in% plus_dual_assigned_qnames)) # store end only assigned reads (insignificant start)
        plus_strand_read_sets[[4]] <- subset(position_data, !(qname %in% combined_assigned_qnames[[1]]) & !(qname %in% combined_assigned_qnames[[2]]) & dna_strand == "+") # store unassigned reads (insignificant start and end)
    }
    else if(i == 4){ # Checking for completion of minus strand read assignment data collection
        minus_qname_identification_data <- rbindlist(combined_position_data, idcol = "Identification_count", fill = TRUE)[qname %in% combined_assigned_qnames[[3]] & qname %in% combined_assigned_qnames[[4]]] # isolate read qnames that have been assigned both a significant start and significant end
        minus_dual_assigned_qnames <- unique(minus_qname_identification_data$qname)
        minus_strand_read_sets[[1]] <- subset(position_data, qname %in% minus_dual_assigned_qnames) # store dual assigned reads (significant start and end)
        minus_strand_read_sets[[2]] <- subset(position_data, qname %in% combined_assigned_qnames[[3]] & !(qname %in% minus_dual_assigned_qnames)) # store start only assigned reads (insignificant end)
        minus_strand_read_sets[[3]] <- subset(position_data, qname %in% combined_assigned_qnames[[4]] & !(qname %in% minus_dual_assigned_qnames)) # store end only assigned reads (insignificant start)
        minus_strand_read_sets[[4]] <- subset(position_data, !(qname %in% combined_assigned_qnames[[3]]) & !(qname %in% combined_assigned_qnames[[4]]) & dna_strand == "-") # store unassigned reads (insignificant start and end)
    }
}

### Isolate the read sets from the sam files, outputting each read set ###
combined_dual_assigned_reads <- rbind(plus_strand_read_sets[[1]], minus_strand_read_sets[[1]])
combined_start_only_assigned_reads <- rbind(plus_strand_read_sets[[2]], minus_strand_read_sets[[2]])
combined_end_only_assigned_reads <- rbind(plus_strand_read_sets[[3]], minus_strand_read_sets[[3]])
combined_unassigned_reads <- rbind(plus_strand_read_sets[[4]], minus_strand_read_sets[[4]])
combined_reads <- list(combined_dual_assigned_reads, combined_start_only_assigned_reads, combined_end_only_assigned_reads, combined_unassigned_reads)
read_set_names <- c("full_length_reads", "5'_end_broken_reads", "3'_end_broken_reads", "unassigned_reads") # dual assigned = full length, start only assigned = 5' end broken, end only assigned = 3' end broken, unassigned = unassigned
for(i in 1:length(combined_reads)){
    for(j in 1:length(sam_files)){
        sam_file <- combined_sam_files[[j]]
        colnames(sam_file) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq")
        sample_reads <- subset(combined_reads[[i]], combined_reads[[i]]$sample == sam_names[j]) # isolate reads for the current sample for the current read set
        sample_qnames <- unique(sample_reads$qname)
        sample_read_set <- subset(sam_file, sam_file$qname %in% sample_qnames) # isolate the same reads from the current sample's sam file
        qual <- rep("*", nrow(sample_read_set))
        sample_read_set <- cbind(sample_read_set, qual) # add a placeholder qual column back to the sam file for output
        colnames(sample_read_set) <- colnames(c(colnames(sample_read_set), "qual"))
        write.table(sample_read_set, file = paste0(sam_names[j], "_", read_set_names[i], ".sam"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) # output the read set as a headerless sam file
    }
}