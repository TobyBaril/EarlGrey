#!/bin/bash


usage()
{
	echo "autoPie.sh -i mergedRepeat.bed -t finalRepeatmasker.tbl -p plotName.pdf -o summaryTable.txt -h show this message"
}

main()
{
	genome_size=$(sed -n '4p' $intbl | rev | cut -f1,1 -d ':' | rev | sed 's/ bp.*//g; s/ //g')
	Rscript ${SCRIPT_DIR}/autoPie.R $inbed ${inbed%.bed}.gff $genome_size $outplot $sumtab
}





## MAIN ##

while getopts i:t:p:o:h option
do
    case "${option}"
        in
        i) inbed=${OPTARG};;
        t) intbl=${OPTARG};;
        p) outplot=${OPTARG};;
	      o) sumtab=${OPTARG};;
        h) usage; exit 1;;
    esac
done

# step 1 - Get genome size and run Rscript

SCRIPT_DIR=INSERTFILEPATHHERE
main
