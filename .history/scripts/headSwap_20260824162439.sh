#!/bin/bash

while getopts i:f:o: option
do
    case "${option}" 
        in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
    esac
done


# 1 - Create Dictionary
SCRIPT_DIR=/data/toby/EarlGrey/scripts

grep ">" "$input" | awk '{print substr($0,2)"\t""ctg_"NR}' > "${input}.dict"
tr -d $'\r' < "${input}.dict" > "${input}.dict.1" && mv "${input}.dict"{.1,}

# 2 - Replace Fasta Headers

${SCRIPT_DIR}/faswap.py ${input}.dict $input > $output
