#!/bin/bash

# program to change filenames of loops to new convention

if [ -z "$1" ]
  then
    echo "must supply input folder"
else
	if ls $1/*_M_*.dat 1> /dev/null 2>&1; then
		echo "step 1, renaming M -> E"
		for f in data/nr/loops/dim_4/K_12/*_M_*.dat; 
			do rename -v 's/_M_/_E_/' $f;
		done
	else
		echo "no files to rename M -> E"
	fi

	if ls $1/*_G_*_B_*.dat 1> /dev/null 2>&1; then
		echo "step 2, combining G and B -> Kappa=G**3*B"
		for f in data/nr/loops/dim_4/K_12/*_G_*_B_*.dat;
			do g=$(echo "$f"|sed -n 's/^.*_G_\([0-9\.]*\)_.*$/\1/p');
			B=$(echo "$f"|sed -n 's/^.*_B_\([0-9\.]*\)_.*$/\1/p');
			kappa=$(python -c "print $g**3*$B");
			fo=$(echo "$f"|sed -n "s/^\(.*\)_G_\([0-9\.]*\)_B_\([0-9\.]\+\)_\(.*\)$/\1_Kappa_$kappa\_\4/p");
			echo $f "renamed as" $fo
			mv $f $fo;
		done
	else
		echo "no files for which to combine G and B -> Kappa=G**3*B"
	fi
fi


