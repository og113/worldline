#!/bin/bash

# program to convert space/tab separated data to a csv

if [ -z "$1" ]
then
	echo "must supply input file";
elif [ -z "$2" ]
then
	outfile=$(echo "$1" | sed 's/\.csv$/.dat/');
else
	outfile=$2;
fi

echo "file in : " $1
echo "file out: " $outfile
#sed -ne 's/"\([0-9\.-]*\)"/\1/gp' $1 >$outfile
cp $1 $outfile
sed -i 's/,/ /g' $outfile
