#!/bin/bash

for i in {4..10..1}
do
	echo "running with K=$i"
	./changeInputs -K $i
	./glmain
	mpirun -n 2 loop
done
