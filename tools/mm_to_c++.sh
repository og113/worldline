#!/bin/bash

wflag=false

# checking if outFile required and getting filename if so
options=':w'
OPTIND=1
while getopts $options option
do
	case $option in
		o  ) wflag=true;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done
shift $((OPTIND-1))

if [ -z "$1" ]
then
	echo "must supply input file"
else
	echo "changing mathematica output to c++ input"
	# generic changes
	sed -i 's/\([1-9]\+\)\*/\1.0\*/g' $1
	sed -i 's/\([1-9]\+\)\.\*/\1.0\*/g' $1
	sed -i 's/\*\([1-9]\+\)/\*\1.0/g' $1
	sed -i 's/\*\([1-9]\+\)\./\*\1.0/g' $1
	sed -i 's/\/\([1-9]\+\)/\/\1.0/g' $1
	sed -i 's/\/\([1-9]\+\)\./\/\1.0/g' $1
	sed -i 's/\([1-9]\+\)\//\1.0\//g' $1
	sed -i 's/\([1-9]\+\)\.\//\1.0\//g' $1
	sed -i 's/Pi/PI/g' $1
	sed -i 's/Power(/pow(/g' $1
	sed -i 's/pow(E,/exp(/g' $1
	sed -i 's/Sqrt(/sqrt(/g' $1
	sed -i 's/\\\[Rho\]/mu/g' $1
	sed -i 's/\\\[Sigma\]/nu/g' $1
	sed -i 's/\\\[Delta\]/d/g' $1
	sed -i 's/\\\[Epsilon\]/a/g' $1
	sed -i 's/\\\[Kappa\]/kappa/g' $1
	sed -i ':a;N;$!ba;s/\(\s*\)\([+-]\)\(\s*\)\n/ \\\n\2 /g' $1
	sed -i 's/^\([+-]\)\(\s*\)/\1 /g' $1
	sed -i 's/^\(\s*\)/ /g' $1
	if $wflag
	then
	# specific changes, for worldline n-r calculations
		sed -i 's/TT(/Dot(/g' $1
		sed -i 's/T(/DX(/g' $1
		sed -i 's/B(\([a-z]\),1 + \([a-z]\),\([a-z]\),1 + \([a-z]\))/B_\1\3/g' $1
		sed -i 's/B(\([a-z]\),-1 + \([a-z]\),\([a-z]\),1 + \([a-z]\))/B_m\1\3/g' $1
		sed -i 's/B(\([a-z]\),1 + \([a-z]\),\([a-z]\),-1 + \([a-z]\))/B_\1m\3/g' $1
		sed -i 's/B(\([a-z]\),-1 + \([a-z]\),\([a-z]\),-1 + \([a-z]\))/B_m\1m\3/g' $1
		sed -i 's/B(\([a-z]*\),\([a-z]*\))/B_\1\2/g' $1
		sed -i 's/B(\([a-z]\),-1 + \([a-z]\))/B_\1m\2/g' $1
		sed -i 's/B(-1 + \([a-z]\),\([a-z]\))/B_m\1\2/g' $1
		sed -i 's/B(-1 + \([a-z]\),-1 + \([a-z]\))/B_m\1m\2/g' $1
		sed -i 's/B(\([a-z]\),1 + \([a-z]\))/B_\1p\2/g' $1
		sed -i 's/B(1 + \([a-z]\),\([a-z]\))/B_p\1\2/g' $1
		sed -i 's/B(1 + \([a-z]\),1 + \([a-z]\))/B_p\1p\2/g' $1
		sed -i 's/EE(\([a-z]*\),\([a-z]*\))/E_\1\2/g' $1
		sed -i 's/EE(\([a-z]\),-1 + \([a-z]\))/E_\1m\2/g' $1
		sed -i 's/EE(-1 + \([a-z]\),\([a-z]\))/E_m\1\2/g' $1
		sed -i 's/EE(-1 + \([a-z]\),-1 + \([a-z]\))/E_m\1m\2/g' $1
		sed -i 's/EE2(\([a-z]*\),\([a-z]*\))/E_\1\2/g' $1
		sed -i 's/EE2(\([a-z]\),-1 + \([a-z]\))/E_\1m\2/g' $1
		sed -i 's/EE2(-1 + \([a-z]\),\([a-z]\))/E_m\1\2/g' $1
		sed -i 's/EE2(-1 + \([a-z]\),-1 + \([a-z]\))/E_m\1m\2/g' $1
		sed -i 's/Dot(\([a-z]\),1 + \([a-z]\),\([a-z]\),1 + \([a-z]\))/T_\1\3/g' $1
		sed -i 's/Dot(\([a-z]\),-1 + \([a-z]\),\([a-z]\),1 + \([a-z]\))/T_m\1\3/g' $1
		sed -i 's/Dot(\([a-z]\),1 + \([a-z]\),\([a-z]\),-1 + \([a-z]\))/T_\1m\3/g' $1
		sed -i 's/Dot(\([a-z]\),-1 + \([a-z]\),\([a-z]\),-1 + \([a-z]\))/T_m\1m\3/g' $1
		sed -i 's/j/y/g' $1
		sed -i 's/k/z/g' $1
		sed -i 's/l/k/g' $1
		sed -i 's/y/i/g' $1
		sed -i 's/z/j/g' $1
		sed -i 's/-1 + \([b-v]\)/m\1/g' $1
		sed -i 's/1 + \([b-v]\)/p\1/g' $1
		sed -i 's/DX(/DX(l,/g' $1
	fi
fi
