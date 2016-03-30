#!/bin/bash
# program to update db/nr/delta.csv with only new results

i='results/nr/nr.csv'
m='db/nr/nr.csv'
o='db/nr/delta.csv'

ls db/nr

comm -23 <(sort $i) <(sort $m) > $o

ls db/nr

i=${i##*/}

ls db/nr

sed -i "s/$/,$i/" $o

ls db/nr
