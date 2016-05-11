#!/bin/bash

wflag=false

# checking if outFile required and getting filename if so
options=':w'
OPTIND=1
while getopts $options option
do
	case $option in
		w  ) wflag=true;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done
shift $((OPTIND-1))

if [ -z "$1" ]
then
	echo "must supply input file"
else
	echo "changing mathematica output to c++ input in $1"
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
	sed -i 's/Cos(/cos(/g' $1
	sed -i 's/Sin(/sin(/g' $1
	sed -i 's/Cosh(/cosh(/g' $1
	sed -i 's/Sinh(/sinh(/g' $1
	sed -i 's/\\\[Rho\]/mu/g' $1
	sed -i 's/\\\[Sigma\]/nu/g' $1
	sed -i 's/\\\[Delta\]/delta/g' $1
	sed -i 's/\\\[Epsilon\]/a/g' $1
	sed -i 's/\\\[Kappa\]/kappa/g' $1
	sed -i 's/\\\[Tau\]/t/g' $1
	sed -i 's/\\\[Beta\]/beta/g' $1
	sed -i ':a;N;$!ba;s/\(\s*\)\([+-]\)\(\s*\)\n/ \\\n\2 /g' $1
	sed -i 's/^\([+-]\)\(\s*\)/\1 /g' $1
	sed -i 's/^\(\s*\)/ /g' $1
	sed -i 's/\/$/\/\\/g' $1
	sed -i 's/\*$/\*\\/g' $1
	# specific changes, for worldline n-r calculations
	if $wflag
	then echo "making specific changes, for worldline n-r calculations"
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
		sed -i 's/Symbolic//g' $1
		sed -i 's/Dd/Dr/g' $1
		sed -i 's/\([,()]\)j/\1yyy/g' $1
		sed -i 's/\([,()]\)j/yyy\1/g' $1
		sed -i 's/\([,()]\)k/\1zzz/g' $1
		sed -i 's/k\([,()]\)/zzz\1/g' $1
		sed -i 's/\([,()]\)l/\1k/g' $1
		sed -i 's/l\([,()]\)/k\1/g' $1
		sed -i 's/yyy/i/g' $1
		sed -i 's/zzz/j/g' $1
		sed -i 's/-1 + \([b-v]\)/m\1/g' $1
		sed -i 's/1 + \([b-v]\)/p\1/g' $1
		sed -i 's/DX(/DX(l,/g' $1
		sed -i 's/Thermak/Thermal/g' $1
		sed -i 's/FThermal(\([a-z]*\),\([a-z]*\))/FThermal_\1\2/g' $1
		sed -i 's/FThermal(\([a-z]\),-1 + \([a-z]\))/FThermal_\1m\2/g' $1
		sed -i 's/FThermal(-1 + \([a-z]\),\([a-z]\))/FThermal_m\1\2/g' $1
		sed -i 's/FThermal(-1 + \([a-z]\),-1 + \([a-z]\))/FThermal_m\1m\2/g' $1
		sed -i 's/Dr(\([a-z]*\),\([a-z]*\))/Dr_\1\2/g' $1
		sed -i 's/Dr(\([a-z]\),-1 + \([a-z]\))/Dr_\1m\2/g' $1
		sed -i 's/Dr(-1 + \([a-z]\),\([a-z]\))/Dr_m\1\2/g' $1
		sed -i 's/Dr(-1 + \([a-z]\),-1 + \([a-z]\))/Dr_m\1m\2/g' $1
		sed -i 's/Dt(\([a-z]*\),\([a-z]*\))/Dt_\1\2/g' $1
		sed -i 's/Dt(\([a-z]\),-1 + \([a-z]\))/Dt_\1m\2/g' $1
		sed -i 's/Dt(-1 + \([a-z]\),\([a-z]\))/Dt_m\1\2/g' $1
		sed -i 's/Dt(-1 + \([a-z]\),-1 + \([a-z]\))/Dt_m\1m\2/g' $1
		sed -i 's/VMeasure(\([a-z]*\),\([a-z]*\))/T_\1\2/g' $1
		sed -i 's/VMeasure(\([a-z]\),-1 + \([a-z]\))/T_\1m\2/g' $1
		sed -i 's/VMeasure(-1 + \([a-z]\),\([a-z]\))/T_m\1\2/g' $1
		sed -i 's/VMeasure(-1 + \([a-z]\),-1 + \([a-z]\))/T_m\1m\2/g' $1
		sed -i 's/d(\([a-z]*\),\([a-z]*\))/r_\1\2/g' $1
		sed -i 's/d(\([a-z]\),-1 + \([a-z]\))/r_\1m\2/g' $1
		sed -i 's/d(-1 + \([a-z]\),\([a-z]\))/r_m\1\2/g' $1
		sed -i 's/d(-1 + \([a-z]\),-1 + \([a-z]\))/r_m\1m\2/g' $1
	fi
fi
