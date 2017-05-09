#!/bin/bash
# program to do a couple of batches of ./shape and ./nrmain runs

for i in {1..20};
	do
	INPUTS="nrinputs/inputs1";
	Epsi=$(python -c "print 0.005*${i}");
	sed -i "s/Epsi\([[:space:]]\+\)\[ \([0-9.]\+\)/Epsi                [ ${Epsi}/" $INPUTS;
	./shape -inputs $INPUTS -shape straightDisjoint -higherOrder 1;
	./nrmain -inputs $INPUTS -pot thermalDisjoint -pass 1 -eigen 1 >> temp/shiftResults.txt;
done;
