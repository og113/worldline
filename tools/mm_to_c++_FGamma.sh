#!/bin/bash

if [ -z "$1" ]
  then
    echo "must supply input file"
else
	echo "changing mathematica output to c++ input, for FGamma"
	sed -i 's/\([1-9]\)\*/\1.0\*/g' $1	
	sed -i 's/\([1-9]\) \([+-]\)/\1.0 \2/g' $1
	sed -i 's/1 + \([a-z]\)/p\1/g' $1
	sed -i 's/\\\[Rho\]/mu/g' $1
	sed -i 's/\\\[Sigma\]/nu/g' $1
	sed -i 's/\\\[Delta\]/delta/g' $1
	sed -i 's/\\\[CapitalGamma\]/G/g' $1
	sed -i 's/Power/pow/g' $1
	sed -i 's/ArcCos/acos/g' $1
	sed -i 's/P(\([a-z]*\))/P[\1]/g' $1
	sed -i 's/Q(\([a-z]*\))/Q[\1]/g' $1
fi
