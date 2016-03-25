#!/bin/bash
# program to simplify results for semiclassical calculation by reducing the number of parameters by one

iflag=false # i for input file
oflag=false # o for output file

options=':i:o:n:'

while getopts $options option
do
	case $option in
		i  ) i=$OPTARG; iflag=true;;
		o  ) o=$OPTARG; oflag=true;;
		n  ) n=$OPTARG;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done

if !($iflag) || !($oflag);
then
    echo "must include options -i -o" >&2
    exit 1
fi

sp='\([[:space:]]\+\)' # a space
spn='\([[:space:]]*\)' # a space, or nor
nm='\([0-9.]\+\)' # a number

if [ $n -eq 0 ];
then
	awk -F "," '{print $1","$2","$3","$4**3*$5","$6","$7","$8","$9","$10","$11}' $i  > $o
else
	grep -vc "^$spn$nm$sp$nm$sp$nm$sp$nm$sp$g$sp$B$sp" $i
	sed -n -e "/^$spn$nm$sp$nm$sp$nm$sp$nm$sp$g$sp$B$sp/!p" $i > $o
fi
