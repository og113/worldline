#!/bin/bash
# program to plot using gp/simple.gp

gpFile='gp/simple.gp'
single=false

if [ -z "$1" ]
then
	echo "must supply input file";
elif [ -z "$2" ]
then
	single=true;
fi

if $single
then
	gnuplot -e "inFile='$*'" $gpFile;
else
	gnuplot -e "inFiles='$*'" $gpFile;
fi
