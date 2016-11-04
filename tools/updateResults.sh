#!/bin/bash
# update the sqlite database with new results
## to do: should get it to take databaseFile and table as inputs (with the given defaults if not input)

if [ -z "$1" ]
then
	echo "must supply input file to updateResults";
	exit 1;
elif [ "$#" -eq 1 ]
then
	echo "must supply table in database";
	exit 1;
	databaseFile="db/nr/nr.sqlite";
elif [ "$#" -eq 2 ]
then
	resultsFile=$1;	
	table=$2;
	databaseFile="db/nr/nr.sqlite";
elif [ "$#" -eq 3 ]
then
	resultsFile=$1;	
	table=$2;
	databaseFile=$3;
else
	echo "too many arguments for updateResults";
	exit 1;
fi

# filenames generated from input
tempResultsFile=$(echo "$resultsFile" | sed 's/results\/nr/temp/');
databaseResultsFile=$(echo "$resultsFile" | sed 's/results/db/');
deltaResultsFile=$(echo "$databaseResultsFile" | sed 's/\([^\/]\+\)$/delta_\1/');

echo "input      : $resultsFile";
echo "updating   : $databaseResultsFile";
echo "           : $deltaResultsFile";
echo "           : $databaseFile - $table";

# copying results into temp file
cat $resultsFile > $tempResultsFile

# getting end part of filename
resultsFileID=${resultsFile##*/}

# adding $resultsFileID to end of temp file
sed -i "s/$/,$resultsFileID/" $tempResultsFile

# putting difference between temp and database files in delta
touch $databaseResultsFile
comm -23 <(sort $tempResultsFile) <(sort $databaseResultsFile) > $deltaResultsFile

deltaLines=$(wc -l < $deltaResultsFile)
echo "new results: $deltaLines";

if [ "$deltaLines" -gt 0 ]
then

# some "here documents" follow. note that there must be no trailing space after the last "HERE"
# this imports results to the sqlite database
sqlite3 $databaseFile <<HERE
.separator ","
.import "$deltaResultsFile" "$table"
HERE

# views to export
view1="s_B_T_a_0.01_K_11_pot_10"
view2="s_B_T_a_0.02_K_11_pot_8"
view3="s_T_a_B_0.01_K_12_pot_10"
view4="s_T_a_B_0.01_K_12_pot_8"
view5="s_sigma_B_T_E_a_0.02_K_11_pot_10"
view6="s_sigma_B_T_E_a_0.02_K_11_pot_8"
echo "updating views ${view1}, ${view2}, ${view3}, ${view4}, ${view5}, ${view6}" 

# exporting views
sqlite3 $databaseFile <<HERE
.mode csv
.separator " "
.header off
.output "temp/${view1}.dat"
SELECT * FROM "$view1";
.output "temp/${view2}.dat"
SELECT * FROM "$view2";
.output "temp/${view3}.dat"
SELECT * FROM "$view3";
.output "temp/${view4}.dat"
SELECT * FROM "$view4";
.output "temp/${view5}.dat"
SELECT * FROM "$view5";
.output "temp/${view6}.dat"
SELECT * FROM "$view6";
HERE

if [ "$?" -eq 0 ]
then
	cat $deltaResultsFile >> $databaseResultsFile
fi


fi


