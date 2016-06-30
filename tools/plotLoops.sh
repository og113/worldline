#!/bin/bash
# program to plot using gp/projection.gp

gpFile='gp/projection.gp'
single=false
cflag=false # should we conevrt to ascii?
oflag=false # is there an outfile?
mflag=false # do we want to print min or max?
bflag=false # do we want to print beta?
kflag=false
K=8

# checking if outFile required and getting filename if so
options=':o:K:mbc'
OPTIND=1
while getopts $options option
do
	case $option in
		o  ) o=$OPTARG; oflag=true;;
		c  ) cflag=true;;
		m  ) mflag=true;;
		b  ) bflag=true;;
		K  ) kflag=true; K=$OPTARG;;
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
	if ! $kflag
	then
		K=$(echo "$@" | sed -n 's/.*K_\([0-9]\+\).*/\1/ p');
		kflag=true;
	fi
	
	if $bflag
	then
		T=$(echo "$@" | sed -n 's/.*T_\([0-9.]\+\).*.dat/\1/ p');
		beta=$(bc <<< "scale=6;1.0/$T");
	fi
	
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
if $cflag
then
	mTemp=""
	for f in $l;
	do fileID=${f##*/};
		tf="data/temp/$fileID";
		./binaryToAscii -b $f -a $tf -loop 1 -K $K
		mTemp+="$tf ";	
	done
	m=${mTemp::-1}
else
	m=$l
fi


# string to plot
if $single
then
	gargs="inFile='$m'"
else
	gargs="inFiles='$m'"	
fi

if $oflag
then
	gargs+="; outFile='$o'"
fi
if $mflag
then
	gargs+="; min='1'; max='1'"
fi
if $bflag
then
	gargs+="; beta=$beta"
fi

# plotting
gnuplot -e "$gargs" $gpFile;
