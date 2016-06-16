#!/bin/bash
# program to do a couple of batches of ./highTemp runs

INPUTS=$(ls nrinputs/highTemp/*)

OUTPUTFOLDER="temp/"

for f in $INPUTS;
	do
	stem=${f##*/};
	fo="${OUTPUTFOLDER}/highTemp_fixDS_0_${stem}";
	fe="${OUTPUTFOLDER}/highTemp_fixDS_0_${stem}_errors";
	./highTemp -inputs $f -fixDS 0 2> $fe 1> $fo;
done;

for f in $INPUTS;
	do
	stem=${f##*/};
	fo="${OUTPUTFOLDER}/highTemp_fixDS_1_${stem}";
	fe="${OUTPUTFOLDER}/highTemp_fixDS_1_${stem}_errors";
	./highTemp -inputs $f -fixDS 1 2> $fe 1> $fo;
done;
