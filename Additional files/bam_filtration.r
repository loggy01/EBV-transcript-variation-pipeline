#!/usr/bin/Rscript --vanilla

######################
# Aaron Logsdon, 2023
######################

# Script that generates n number of filtered bam files from n number of input bam files
# Example use is to filter bam files for reads aligning to EBV, retaining only reads that align to EBV (Rscript ./bam_filtration.r ./input.bam name_one 3 ./EBV_reference_name.txt TRUE)

library(Rsamtools)

args <- commandArgs(trailingOnly = TRUE)
bam_files <- strsplit(args[1], ",")[[1]] # user bam files (list of directories)
bam_names <- strsplit(args[2], ",")[[1]] # names matching the bam files (list of strings)
data_index <- as.numeric(args[3]) # index representing the column in each bam file that will be filtered against (see https://samtools.github.io/hts-specs/SAMv1.pdf) (integer)
data_charateristics <- strsplit(args[4], ",")[[1]] # single column (with column header) text files matching the bam files. The data here is applied across the specified data_index for the matching bam file (list of directories)
comparison <- args[5] # TRUE or FALSE. If TRUE, the bam files will be filtered to include only reads that match the data in the data_charateristics file. If FALSE, the bam files will be filtered to include only reads that do not match the data in the data_charateristics file (boolean)

for(i in 1:length(bam_files)){
    input_bam <- BamFile(bam_files[i])
    output_bam <- paste0(bam_names[i], "_filtered.bam")
    read_filter_rules <- ""
    current_data_charateristics <- read.table(data_charateristics[i], header = TRUE, sep = "\t")
    current_data_charateristics <- current_data_charateristics[, 1]
    if(comparison == "TRUE"){
        read_filter_rules <- FilterRules(list(keepRead = function(read) read[[data_index]] %in% current_data_charateristics == TRUE))
    }
    else if(comparison == "FALSE"){
        read_filter_rules <- FilterRules(list(keepRead = function(read) read[[data_index]] %in% current_data_charateristics != TRUE))
    }
    filterBam(input_bam, output_bam, indexDestination = TRUE, filter = read_filter_rules)
}
