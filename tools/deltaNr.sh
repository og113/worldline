#!/bin/bash
# program to update db/nr/delta.csv with only new results

i='results/nr/nr.csv'
m='db/nr/nr.csv'
o='db/nr/delta.csv'

comm -23 <(sort $i) <(sort $m) > $o

i=${i##*/}

sed -ie "s/$/,$i/" $o
