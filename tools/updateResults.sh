#!/bin/bash
# update the sqlite database with new results

if [ -z "$1" ]
then
	echo "must supply input file to updateResults";
	exit 1;
elif [ "$#" -eq 1 ]
then
	echo "must supply table in database";
	exit 1;
	databaseFile="db/nr.sqlite";
elif [ "$#" -eq 2 ]
then
	resultsFile=$1;	
	table=$2;
	databaseFile="db/nr.sqlite";
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
tempResultsFile=$(echo "$resultsFile" | sed 's/results/temp/');
databaseResultsFile=$(echo "$resultsFile" | sed 's/data\/results/db/');
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

if [ "$?" -eq 0 ]
then
	cat $deltaResultsFile >> $databaseResultsFile
fi

fi

# views to export
views=""
if [ "$table" == "nr6" ];
then
	views="s_sigma_B_T_a_E_K_11_pot_13_nr6 s_sigma_B_T_a_E_K_11_pot_15_nr6 kta_pot_13 kta_pot_15"
	#"s_B_T_a_0.01_K_11_pot_10 s_B_T_a_0.02_K_11_pot_8 s_T_a_B_0.01_K_12_pot_10 s_T_a_B_0.01_K_12_pot_8 s_sigma_B_T_E_a_0.02_K_11_pot_10 s_sigma_B_T_E_a_0.02_K_11_pot_8 s_sigma_B_T_a_E_K_11_pot_8 s_sigma_B_T_a_E_K_11_pot_10 s_sigma_B_T_a_E_K_11_pot_13 s_sigma_B_T_a_E_K_11_pot_14 s_sigma_B_T_a_E_K_11_pot_15 kta_pot_8 kta_pot_10 kta_pot_12 kta_pot_13 kta_pot_14 kta_pot_15"
elif [ "$table" == "nr7" ];
then
	views="s_sigma_B_T_a_E_K_12_pot_13_nr7 s_sigma_B_T_a_E_K_12_pot_15_nr7 s_sigma_B_T_a_E_K_12_pot_16_nr7 kta_pot_13_K_12 kta_pot_15_K_12 s_T_a_B_0.1_K_12_pot_13 s_T_a_B_0.1_K_12_pot_15  s_T_a_B_0.5_K_12_pot_13 s_T_a_B_0.5_K_12_pot_15 s_B_a_K_12_pot_1 s_B_a_K_12_pot_0"
elif [ "$table" == "eigs7" ];
then
	views="eigs7straight_pot_15"
elif [ "$table" == "eigs7temp" ];
then
	views="eigs7tempstraight_pot_15"
elif [ "$table" == "eigs8" ];
then
	views="eigs7tempstraight_pot_15"
fi
	
echo "updating views:"

# exporting views
for v in $views;
do
echo $v
sqlite3 $databaseFile <<HERE
.mode csv
.separator " "
.header off
.output "db/nr/plotsHand/${v}.dat"
SELECT * FROM "${v}";
HERE
done

