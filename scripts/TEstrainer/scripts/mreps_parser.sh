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
exit 1
fi

set -eo pipefail
if mreps -minsize 20 -maxsize 50000 -minperiod 1 -maxperiod 5000 -exp 2 -fasta "$INFILE" > "${INFILE}.mreps.raw" 2> "${INFILE}.mreps.stderr"; then
    cat "${INFILE}.mreps.stderr" >&2
else
    status=$?
    cat "${INFILE}.mreps.stderr" >&2
    if [ "$status" -eq 17 ] && grep -Fxq "Error: Too many N's in the window" "${INFILE}.mreps.stderr"; then
        echo "WARNING: mreps did not finish for ${INFILE} due to ambiguous bases; its repeat results may be incomplete. TRF and SA-SSR results are still used." >&2
    else
        exit "$status"
    fi
fi
awk '/->/ && !/from/' "${INFILE}.mreps.raw" | sed 's/ \+/\t/g' | sed 's/\t\+/\t/g' | sed 's/\t//10g' | cut -f2,4,6,7,8,9,10 | sed 's/ ://g' | sed 's/>//g' | sed 's/<//g' | sed 's/ //g' | sed 's/\[//g' | sed 's/\]//g' > ${INFILE}.mreps.temp
awk '{OFS="\t"}{print FILENAME,$1,$2,$3,$4,$5,$6,$7}' ${INFILE}.mreps.temp | sed 's/.*\///;s/.fasta.mreps.temp//' > ${INFILE}.mreps
rm "${INFILE}.mreps.temp" "${INFILE}.mreps.raw" "${INFILE}.mreps.stderr"
