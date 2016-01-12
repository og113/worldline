#!/bin/bash

if [ -z "$1" ]
  then
    echo "must supply input file"
else
	echo "changing mathematica output to c++ input"
	sed -i 's/\\\[Rho\]/m/g' $1
	sed -i 's/\\\[Sigma\]/n/g' $1
	sed -i 's/\\\[Delta\]/d/g' $1
	sed -i 's/TT(/Dot(/g' $1
	sed -i 's/T(/DX(/g' $1
	sed -i 's/Power(/pow(/g' $1
	sed -i ':a;N;$!ba;s/\(\s*\)\([+-]\)\(\s*\)\n/ \\ \n\2 /g' $1
	sed -i 's/^\([+-]\)\(\s*\)/\1 /g' $1
	sed -i 's/^\(\s*\)/ /g' $1
	sed -i 's/B(\([a-z]\),1 + \([a-z]\),\([a-z]\),1 + \([a-z]\))/B_\1\3/g' $1
	sed -i 's/B(\([a-z]\),-1 + \([a-z]\),\([a-z]\),1 + \([a-z]\))/B_m\1\3/g' $1
	sed -i 's/B(\([a-z]\),1 + \([a-z]\),\([a-z]\),-1 + \([a-z]\))/B_\1m\3/g' $1
	sed -i 's/B(\([a-z]\),-1 + \([a-z]\),\([a-z]\),-1 + \([a-z]\))/B_m\1m\3/g' $1
fi
