#!/bin/bash

results="results/nr/nr6.csv results/nr/nr6error.csv results/nr/nr7.csv results/nr/nr7error.csv"

echo "removing lines with <<<, >>> and === from:"

for f in $results;
do
	echo $f
	sed -i '/<<</d' $f
	sed -i '/>>>/d' $f
	sed -i '/===/d' $f
done
