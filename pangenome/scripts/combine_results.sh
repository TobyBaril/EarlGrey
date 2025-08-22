#!/bin/bash

# Usage: combine_results.sh file1 file2 ... fileN output_file

# Get the output file (last argument)
output_file="${@: -1}"

# Get all input files (all but last argument)
input_files=("${@:1:$#-1}")

# Clear the output file if it exists
> "$output_file"

# Combine the results from the input files
cat "${input_files[@]}" >> "$output_file"