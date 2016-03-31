#!/bin/bash
# program to update db/nr/delta.csv with only new results

raw='results/nr/nr.csv'
rawf='temp/nrf.csv'
db='db/nr/nr.csv'
delta='db/nr/delta.csv'

cat $raw > $rawf

raw=${raw##*/}

sed -i "s/$/,$raw/" $rawf

comm -23 <(sort $rawf) <(sort $db) > $delta

# to append db/nr/delta.csv to the end of db/nr/nr.csv, type "cat db/nr/delta.csv >> db/nr/nr.csv"
