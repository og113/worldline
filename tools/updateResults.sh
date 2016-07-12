#!/bin/bash
# update the sqlite database with new results

if [ -z "$1" ]
then
	echo "must supply input file to updateResults";
	exit 1;
elif [ -z "$2" ]
then
	databaseFile="db/nr/nr.sqlite"
	table="nrtemp"
	
	resultsFile=$1;
	tempResultsFile=$(echo "$1" | sed 's/results\/nr/temp/');
	databaseResultsFile=$(echo "$1" | sed 's/results/db/');
	deltaResultsFile=$(echo "$databaseResultsFile" | sed 's/\([^\/]\+\)$/delta_\1/');
	
	echo "input     : $resultsFile";
	echo "updating  : $databaseResultsFile";
	echo "          : $deltaResultsFile";
	echo "          : $databaseFile";
	
else
	echo "too many arguments for updateResults";
	exit 1;
fi

# copying results into temp file
cat $resultsFile > $tempResultsFile

# getting end part of filename
resultsFileID=${resultsFile##*/}

# adding $resultsFileID to end of temp file
sed -i "s/$/,$resultsFileID/" $tempResultsFile

# putting difference between temp and database files in delta
touch $databaseResultsFile
comm -23 <(sort $tempResultsFile) <(sort $databaseResultsFile) > $deltaResultsFile

# a "here document" follows. note that there must be no trailing space after the last "HERE"
# this iports results to the sqlite database
sqlite3 $databaseFile <<HERE
.separator ","
.import "$deltaResultsFile" "$table"
HERE

cat $deltaResultsFile >> $databaseResultsFile
