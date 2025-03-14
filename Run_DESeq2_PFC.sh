#!/bin/bash

region="PFC"

# Pseudobulk file location
input_dir="/home/AD/a5green/scorch/DESeq2/pb_files_2/"

# Date on pseudobulk file
pb_date="022525"

output_dir="/home/AD/a5green/scorch/DESeq2/deseq_results"
# Get current date in YYYYMMDD format
current_date=$(date +"%Y%m%d")
# Append date to output directory
output_dir_2="${output_dir}_${region}_${current_date}"
# Create output directory if it doesn't exist
mkdir -p "$output_dir_2"

# Run the R script with arguments
Rscript DESeq2_Script.R "$region" "$input_dir" "$pb_date" "$output_dir_2"
