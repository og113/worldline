#!/bin/bash

# program to change filenames of loops to new convention

if [ -z "$1" ]
  then
    echo "must supply input folder"
else
	test=$(ls $1/*_Kappa_*.dat | wc -l)
	if (("$test" > "0")); then
		echo "step 1, renaming Kappa -> kappa"
		for f in $1/*_Kappa_*.dat;
			do fo=$(echo "$f"|sed "s/_Kappa_/_kappa_/")
			echo $f "renamed as" $fo
			mv $f $fo;
		done
	else
		echo "no files to rename Kappa -> kappa"
	fi
fi
