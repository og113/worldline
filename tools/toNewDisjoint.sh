#!/bin/bash

if [ -z "$1" ]
then
	echo "must supply input file"
else
	echo "changing loop function to disjoint function in $1"
	sed -i 's/posNeigh(/posNeighDisjoint(/g' $1
	sed -i 's/negNeigh(/negNeighDisjoint(/g' $1
	sed -i 's/Dot(\([^)]\+\))/DotDisjoint(\1,beta)/g' $1
	sed -i 's/DX(\([^,]\+\),\([^,]\+\),\([^,)]\+\))/DXDisjoint(\1,\2,\3,beta)/g' $1
fi
