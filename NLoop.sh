#!/bin/bash

for i in {0..1000..200}
do
	echo "running with Nms=$i"
	./changeInputs -Nms $i
	./glmain
	mpirun -n 2 loop
done
