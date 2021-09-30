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

grep ">" $input | awk '{printf("ctg_%01d %s\n", NR, $0)}'  | sed 's/ /\t/g; s/>//g' | awk '{OFS="\t"}{print $2, $1}' > ${input}.dict
tr -d $'\r' < ${input}.dict > ${input}.dict.1 && mv ${input}.dict{.1,}

# 2 - Replace Fasta Headers

faswap.py ${input}.dict $input > $output
