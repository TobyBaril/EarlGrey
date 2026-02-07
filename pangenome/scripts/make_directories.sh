#!/bin/bash

directory=$1

mkdir -p ${directory}/${species}_EarlGrey/ && OUTDIR=${directory}/${species}_EarlGrey
if [ ! -z "$RepSpec" ] || [ ! -z "$startCust" ] ; then
	mkdir -p $OUTDIR/${species}_RepeatMasker/
fi
mkdir -p $OUTDIR/${species}_Database/
mkdir -p $OUTDIR/${species}_RepeatModeler/
mkdir -p $OUTDIR/${species}_strainer/
mkdir -p $OUTDIR/${species}_RepeatMasker_Against_Custom_Library/
mkdir -p $OUTDIR/${species}_RepeatLandscape/
mkdir -p $OUTDIR/${species}_mergedRepeats/
mkdir -p ${OUTDIR}/${species}_summaryFiles/
if [ ! -z "$heli" ]; then
	mkdir -p ${OUTDIR}/${species}_heliano/
fi
