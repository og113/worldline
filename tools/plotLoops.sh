#!/bin/bash
# program to plot using gp/projection.gp

gpFile='gp/projection.gp'
single=false

if [ -z "$1" ]
then
	echo "must supply input file";
elif [ -z "$2" ]
then
	single=true;
fi

m=""
for f in "$@";
	do fileID=${f##*/};
	tf="data/temp/$fileID";
	./binaryToAscii -b $f -a $tf -loop 1
	m+="$tf";
	if !($single)
	then
		m+=" "
	fi
done

if $single
then
	gnuplot -e "inFile='$m'" $gpFile;
else
	gnuplot -e "inFiles='$m'" $gpFile;
fi
