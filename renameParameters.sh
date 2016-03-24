#!/bin/bash

# program to change filenames of loops to new convention
for f in data/nr/loops/dim_4/K_12/*.dat; 
	do rename -v 's/_M_/_E_/' $f;
done

for f in data/nr/loops/dim_4/K_12/*.dat;
do
g=$(echo "$f"|sed -n 's/^.*_G_\([0-9\.]*\)_.*$/\1/p');
B=$(echo "$f"|sed -n 's/^.*_B_\([0-9\.]*\)_.*$/\1/p');
kappa=$(python -c "print $g**3*$B");
fo=$(echo "$f"|sed -n "s/^\(.*\)_G_\([0-9\.]*\)_B_\([0-9\.]*\)_\(.*\)$/\1_Kappa_$kappa\_\4/p");
echo $f "renamed as" $fo
mv $f $fo;
done
