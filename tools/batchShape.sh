#!/bin/bash
# program to do a couple of batches of ./shape and ./nrmain runs

if [ $# -ne 1 ]; then
    INPUTS="nrinputs/inputs0";
else
	INPUTS=$1;
fi

for i in {1..2};
	do
#	Epsi=$(python -c "print 0.005*${i}");
	Epsi=$(bc <<< "scale=3;0.10*${i}");
	echo "a = 0${Epsi}";
	sed -i "s/Epsi\([[:space:]]\+\)\[ \([0-9.]\+\)/Epsi                [ 0${Epsi}/" $INPUTS;
	./shape -base /slow/data/dc-goul1/worldline/ -inputs $INPUTS -shape straightDisjoint -higherOrder 0;
	./nrmain -base /slow/data/dc-goul1/worldline/ -inputs $INPUTS -guess 1 -step 0 -old 0 -pot thermalDisjoint -gaussian 0 -gaussianLR 0 -mu_a 0 -nonrel 0 -allTests 0 -fixdz 1 -fixtlr 1 -fixdislr 1 -print none -auto_a 0 -conservation 0 -redo 1 -redoErrors 0 -verbose 1 -pass 1 -eigen 1
	# N.B. should have -pot thermalDisjointLR, but just doing a test run without
done;
