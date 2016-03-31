#!/bin/bash
# program to plot using gp/projection.gp

gpFile='gp/projection.gp'
single=false
oflag=false

# checking if outFile required and getting filename if so
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

# checking that input files were specified
if [ -z "$1" ]
	then
	echo "must supply input file";
fi

lTemp=""
for f in "$@";
do 
	if [ -e "$f" ] && [ -f "$f" ];
	then
		lTemp+="$f ";
	elif [ -e "$f" ] && [ -d "$f" ]; # only looks one directory deep, and fails if there are subdirectories
	then
		lTemp+=$(ls $f/*.dat | tr '\n' ' ')
	fi
done
l=${lTemp::-1}

# checking if multiple files are being plotted
lwords=$(echo $l | wc -w)
if [ "$lwords" -gt 1 ]
then
	single=false;
fi

# converting binary to ascii files
mTemp=""
for f in $l;
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
