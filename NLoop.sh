#!/bin/bash

for i in {0..24..1}
do
	g=$(awk "BEGIN {printf \"%.2f\",$i/4}")
	echo "running with Nms = $g"
	./changeInputs -g $i
	./glmain
	mpirun -n 2 loop
done
