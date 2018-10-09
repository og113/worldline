#!/bin/bash

iflag=false # i for input file
oflag=false # o for output file

options=':i:o:'
OPTIND=1
while getopts $options option
do
	case $option in
		i  ) i=$OPTARG; iflag=true;;
		o  ) o=$OPTARG; oflag=true;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done
shift $((OPTIND-1))

if !($iflag) || !($oflag);
then
    echo "must include options -i -o" >&2
    exit 1
fi
