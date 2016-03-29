#!/bin/bash
# program to plot using gp/projection.gp

gpFile='gp/projection.gp'
single=false
oflag=false

# checking of outFile required and getting filename if so
options=':o:'
OPTIND=1
while getopts $options option
do
	case $option in
		o  ) o=$OPTARG; oflag=true;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done
shift $((OPTIND-1))

# checking if multiple files are being plotted
if [ -z "$1" ]
then
	echo "must supply input file";
elif [ -z "$2" ]
then
	single=true;
fi

# converting binary to ascii files
mTemp=""
for f in "$@";
	do fileID=${f##*/};
	tf="data/temp/$fileID";
	./binaryToAscii -b $f -a $tf -loop 1
	mTemp+="$tf ";
done
m=${mTemp::-1}

# plotting
if $single
then
	if $oflag
	then
		gnuplot -e "inFile='$m'; outFile='$o'" $gpFile;
	else
		gnuplot -e "inFile='$m'" $gpFile;
	fi
else
	if $oflag
	then
		gnuplot -e "inFiles='$m'; outFile='$o'" $gpFile;
	else
		gnuplot -e "inFiles='$m'" $gpFile;
	fi
fi
