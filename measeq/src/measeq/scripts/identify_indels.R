#!/usr/bin/env Rscript
# R code that identifies INDELS
# Load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
# Retrieve arguments from bash
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
results <- args[2]
# Set working directory
setwd(wd)
# Specify path to unzipped vcf folders
folder_path <- paste0(results, "/variants/ivar")
# Create a list of all the vcf files in the folder to be read 
vcf_files <- list.files(folder_path, pattern = "\\.vcf$", full.names = TRUE)
# Create an empty list for the data to be appended to
filtered_data <- list()
# Run a for loop reading the vcf files and filtering based on the number of characters in the REF or ALT columns 
for (file in vcf_files){
    data <- fread(file, skip=15)
    data$name <- file
    filtered_rows <- data[nchar(data$REF) > 1 | nchar(data$ALT) > 1,]
    filtered_data[[file]] <- filtered_rows
}
# Filter out vcf files that don't have indels
filtered_list <- Filter(function(df) nrow(df) > 0, filtered_data)
# Combine all the files with indels into one dataframe 
combined_indels <- rbindlist(filtered_list, use.names = FALSE)
# Remove the full path in the name column so that the sample for each indel can easily be identified 
combined_indels$name <- gsub(folder_path, "", combined_indels$name)
# Only keep rows where there are additions or deletions
combined_indels <- subset(combined_indels, nchar(REF) !=nchar(ALT))
combined_indels <- combined_indels[, c(1:2, 4:5, 11)]
# Write file 
write.csv(combined_indels, file = paste0(folder_path,"/indels.csv"), row.names = FALSE)
