#!/bin/bash

if [ -z "$1" ]
  then
    echo "must supply input file"
else
	outfile=$(echo "$1" | sed 's/\.\([a-z]*\)$/.csv/')
	echo "file in : " $1
	echo "file out: " $outfile
	sed -ne 's/\([0-9\.]\+\)\([[:space:]]\+\)/\1,/gp' $1>$outfile
fi
