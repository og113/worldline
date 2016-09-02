#!/bin/bash

if [ -z "$1" ]
  then
    echo "must supply input file"
else
	echo "changing mathematica output to c++ input, for FGamma"
	sed -i 's/\([1-9]\+\)\*/\1.0\*/g' $1
	sed -i 's/\([1-9]\+\)\.\*/\1.0\*/g' $1
	sed -i 's/\*\([1-9]\+\)/\*\1.0/g' $1
	sed -i 's/\*\([1-9]\+\)\./\*\1.0/g' $1
	sed -i 's/\/\([1-9]\+\)/\/\1.0/g' $1
	sed -i 's/\/\([1-9]\+\)\./\/\1.0/g' $1
	sed -i 's/\([1-9]\+\)\//\1.0\//g' $1
	sed -i 's/\([1-9]\+\)\.\//\1.0\//g' $1
	sed -i 's/1 + \([a-z]\)/p\1/g' $1
	sed -i 's/Pi/PI/g' $1
	sed -i 's/\\\[Beta\]/beta/g' $1
	sed -i 's/\\\[Rho\]/mu/g' $1
	sed -i 's/\\\[Sigma\]/nu/g' $1
	sed -i 's/\\\[Delta\]/delta/g' $1
	sed -i 's/\\\[CapitalGamma\]/G/g' $1
	sed -i 's/Power/pow/g' $1
	sed -i 's/pow(E,/exp(/g' $1
	sed -i 's/Sqrt(/sqrt(/g' $1
	sed -i 's/ArcCos/acos/g' $1
	sed -i 's/Cos(/cos(/g' $1
	sed -i 's/Sin(/sin(/g' $1
	sed -i 's/Cosh(/cosh(/g' $1
	sed -i 's/Sinh(/sinh(/g' $1
	sed -i 's/EllipticTheta(3,(PI\*t)/jbdtheta3(t/g' $1
	sed -i 's/EllipticThetaPrime(3,(PI\*t)/jbderivutheta3(t/g' $1
	sed -i 's/Derivative(0,1,0)(EllipticThetaPrime)(3,(PI\*t)/jbderiv2utheta3(t/g' $1
	sed -i 's/P(\([a-z]*\))/P[\1]/g' $1
	sed -i 's/Q(\([a-z]*\))/Q[\1]/g' $1
fi
