#!/bin/bash
# program to plot using gp/simple.gp

gpFile='gp/simple.gp'
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

# plotting
if $single
then
	if $oflag
	then
		gnuplot -e "inFile='$*'; outFile='$o'" $gpFile;
	else
		gnuplot -e "inFile='$*'" $gpFile;
	fi
else
	if $oflag
	then
		gnuplot -e "inFiles='$*'; outFile='$o'" $gpFile;
	else
		gnuplot -e "inFiles='$*'" $gpFile;
	fi
fi
