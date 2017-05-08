#!/bin/bash
# program to do a couple of batches of ./shape and ./nrmain runs

for i in {0..10};
	do
	INPUTS="nrinputs/shape/inputs${i}";
#	Mu=$(python -c "print 0.1*${i}");
	./shape -inputs $INPUTS -shape straightDisjoint -higherOrder 1 -shift 1;
	./nrmain -inputs $INPUTS -pot thermalDisjoint -pass 1 -eigen 1 >> temp/shiftResults.txt;
done;
