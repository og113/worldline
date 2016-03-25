#!/bin/bash
# program to simplify results for semiclassical calculation by reducing the number of parameters by one

gflag=false
Bflag=false
iflag=false # i for input file
oflag=false # o for output file

options=':g:B:i:o:n:'

while getopts $options option
do
	case $option in
		g  ) g=$OPTARG; gflag=true;;
		B  ) B=$OPTARG; Bflag=true;;
		i  ) i=$OPTARG; iflag=true;;
		o  ) o=$OPTARG; oflag=true;;
		n  ) n=$OPTARG;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done

if !($gflag) || !($Bflag) || !($iflag) || !($oflag);
then
    echo "must include options -g -B -i -o" >&2
    exit 1
fi

if [ $n -eq 0 ];
then
	grep -vc "^\([[:space:]]*\)\([0-9.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)$g\([[:space:]]\+\)$B\([[:space:]]\+\)" $i
	sed -n -e "/^\([[:space:]]*\)\([0-9.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)$g\([[:space:]]\+\)$B\([[:space:]]\+\)/!p" $i > $o
else
	grep -vc "^\([[:space:]]*\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)$g\([[:space:]]\+\)$B\([[:space:]]\+\)" $i
	sed -n -e "/^\([[:space:]]*\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)\([0-9\.]\+\)\([[:space:]]\+\)$g\([[:space:]]\+\)$B\([[:space:]]\+\)/!p" $i > $o
fi
