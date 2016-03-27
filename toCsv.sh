#!/bin/bash

# program to convert space/tab separated data to a csv

if [ -z "$1" ]
then
	echo "must supply input file";
elif [ -z "$2" ]
then
	outfile=$(echo "$1" | sed 's/\.\([a-z]*\)$/.csv/');
else
	outfile=$2;
fi

echo "file in : " $1
echo "file out: " $outfile
sed -ne 's/^\([[:space:]]*\)\([0-9\.]\+\)/\2/gp' $1>$outfile
sed -ni 's/\([0-9\.]\+\)\([[:space:]]\+\)/\1,/gp' $outfile
