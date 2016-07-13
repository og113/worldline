#!/bin/bash
# update the sqlite database with new results
## to do: should get it to take databaseFile and table as inputs (with the given defaults if not input)

if [ -z "$1" ]
then
	echo "must supply input file to updateResults";
	exit 1;
elif [ "$#" -eq 1 ]
then
	resultsFile=$1;	
	table="nrtemp";
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

echo "input     : $resultsFile";
echo "updating  : $databaseResultsFile";
echo "          : $deltaResultsFile";
echo "          : $databaseFile";

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

if [ "$deltaLines" -gt 0 ]
then

# a "here document" follows. note that there must be no trailing space after the last "HERE"
# this iports results to the sqlite database
sqlite3 $databaseFile <<HERE
.separator ","
.import "$deltaResultsFile" "$table"
HERE

cat $deltaResultsFile >> $databaseResultsFile

fi


