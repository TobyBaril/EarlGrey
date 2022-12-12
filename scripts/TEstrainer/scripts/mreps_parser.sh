#!/bin/bash
# run mreps and parse to usable format

usage() { echo "Usage: [-i input sequence (fasta file)] [-h Print this help]" 1>&2; exit 1; }

while getopts i:h flag
do
    case "${flag}" in
        i) INFILE=${OPTARG};;
        h | *)
          usage
          exit
    esac
done

if [ ! -f "${INFILE}" ]; then
echo "Error: file not found"
exit
fi

mreps -minsize 20 -maxsize 50000 -minperiod 1 -maxperiod 5000 -exp 2 -fasta $INFILE | grep "\->" | grep -v "from" | sed 's/ \+/\t/g' | sed 's/\t\+/\t/g' | sed 's/\t//10g' | cut -f2,4,6,7,8,9,10 | sed 's/ ://g' | sed 's/>//g' | sed 's/<//g' | sed 's/ //g' | sed 's/\[//g' | sed 's/\]//g' > ${INFILE}.mreps.temp
awk '{OFS="\t"}{print FILENAME,$1,$2,$3,$4,$5,$6,$7}' ${INFILE}.mreps.temp | sed 's/.*\///;s/.fasta.mreps.temp//' > ${INFILE}.mreps
rm ${INFILE}.mreps.temp
