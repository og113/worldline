#!/bin/bash
# program to update results/nr/delta.csv with only new results

i='results/nr/nr.csv'
m='db/nr/nr.csv'
o='db/nr/delta.csv'

comm -23 <(sort $i) <(sort $m) > $o

i=${i##*/}

sed -i "s/$/,$i/" $o

# to append db/nr/delta.csv to the end of db/nr/nr.csv, type "cat db/nr/delta.csv >> db/nr/nr.csv"

ip='results/nr/nr_pass.csv'
mp='db/nr/nr_pass.csv'
op='db/nr/delta_pass.csv'

comm -23 <(sort $ip) <(sort $mp) > $op

ip=${i##*/}

sed -i "s/$/,$ip/" $op
