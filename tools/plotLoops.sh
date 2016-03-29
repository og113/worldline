#!/bin/bash
# program to plot using gp/projection.gp

gpFile='gp/projection.gp'
many=false

if [ -z "$1" ]
then
	echo "must supply input file";
elif [ -z "$2" ]
then
	many=true;
fi

if $many
then
	gnuplot -e "inFile='$1'" $gpFile;
else
	gnuplot -e "inFiles='$*'" $gpFile;
fi
