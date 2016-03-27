#!/bin/bash
# program to simplify results for semiclassical calculation by reducing the number of parameters by one

iflag=false # i for input file
oflag=false # o for output file

options=':i:o:'

while getopts $options option
do
	case $option in
		i  ) i=$OPTARG; iflag=true;;
		o  ) o=$OPTARG; oflag=true;;
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

# for nrmain_cosmos_8.csv
#awk -F "," '{print $1","$2","$3","$4","$5**3*$6","$7","$8","$9","$10","$11","$12","$13",,,,,,"$14","$15","$16","$17","$18","$19","$20","$21}' $i  > $o

# for nrmain_laptop_7.csv and nrmain_cosmos_7.dat
#awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$10","$11",,,,,,"$12","$13","$14","$15","$16","$17","$18","$19}' $i  > $o

# for nrmain_cosmos_6.dat
#awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$10","$11",,,,,,"$12","$13","$14","$15","$16","$17","$18","}' $i  > $o

# for nrmain_cosmos_5.dat
#awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$10","$11",,,,,,"$12","$13","$14","$15","$16","$17",,"}' $i  > $o

# for nrmain_cosmos_4.dat
#awk -F "," '{print $1","$2",0,"$3","$4**3*$5",,"$6","$7",,"$8","$9","$10",,,,,,"$11","$12","$13","$14","$15","$16",,"}' $i  > $o

# for nrmain_cosmos_3.dat
#awk -F "," '{print $1","$2",0,"$3","$4**3*$5",,"$6",1.0,,"$7","$8",,,,,,,"$9","$10",,,"$11","$12",,"}' $i  > $o

# for nrmain_cosmos_2.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5",,"$6",1.0,,"$7","$8",,,,,,,"$9","$10",,,"$11",,,"}' $i  > $o

# for nrmain_cosmos.dat
#awk -F "," '{print $1","$2",0,"$3","$4**3*$5",0,"$6",1,,"$7","$10",,"$8","$9",,,,"$11",,,,,,,"}' $i  > $o
