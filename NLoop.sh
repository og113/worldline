#!/bin/bash

for i in {0..30..1}
do
	g=$(awk "BEGIN {printf \"%.2f\",$i/6}")
	echo "running with g = $g"
	./changeInputs -g $g
	./glmain
	mpirun -n 2 loop
done
