#!/bin/bash

#counting occurances of ',' on first few lines of file

lines=10

if [ -z "$1" ]
then
	echo "must supply input file";
else
	awk -F ','  'NR<='$lines'{print NF-1, NR}' $1
fi


