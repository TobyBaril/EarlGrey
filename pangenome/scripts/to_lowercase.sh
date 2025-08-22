#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

input_file=$1
output_file=$2

# Convert the contents of the input file to lowercase and save to the output file
tr '[:upper:]' '[:lower:]' < "$input_file" > "$output_file"