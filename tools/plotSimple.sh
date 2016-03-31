#!/bin/bash
# program to plot using gp/simple.gp

gpFile='gp/simple.gp'
single=true
oflag=false
tflag=false

# checking if outFile required and getting filename if so
options=':o:t:'
OPTIND=1
while getopts $options option
do
	case $option in
		o  ) o=$OPTARG; oflag=true;;
		o  ) t=$OPTARG; tflag=true;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done
shift $((OPTIND-1))

# checking that input files were specified
if [ -z "$1" ]
	then
	echo "must supply input file";
fi

mTemp=""
for f in "$@";
do 
	if [ -e "$f" ] && [ -f "$f" ];
	then
		mTemp+="$f ";
	elif [ -e "$f" ] && [ -d "$f" ]; # only looks one directory deep, and fails if there are subdirectories
	then
		tflag=true
		t=$f
		mTemp+=$(ls $f/*.dat | tr '\n' ' ')
		mTemp+=" "
	fi
done
m=${mTemp::-1}

# checking if multiple files are being plotted
mwords=$(echo $m | wc -w)
if [ "$mwords" -gt 1 ]
then
	single=false;
fi

# plotting
if $single
then
	gargs="inFile='$m'"
	if $oflag
	then
		gargs+="; outFile='$o'"
	fi
	if $tflag
	then
		gargs+="; Title='$t'"
	fi
	
	gnuplot -e "$gargs" $gpFile;

else
	gargs="inFiles='$m'"
	if $oflag
	then
		gargs+="; outFile='$o'"
	fi
	if $tflag
	then
		gargs+="; Title='$t'"
	fi
	
	gnuplot -e "$gargs" $gpFile;
	
fi
