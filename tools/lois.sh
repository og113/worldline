#!/bin/bash

# a quick and dirty bash program to increase "B" and submit a job "subs/nrmain0.sub" with the appropriate arguments

# defining some files
resultsFile="results/nr/nr4.csv"
inputsFile="nrinputs/highTemp/inputs0"
#inputsFile=inputs4
jobFile="subs/nrmain0.sub"

echo "results file = $resultsFile"
echo "inputs file = $inputsFile"
echo "job file = $jobFile"

# finding the parameter "B" on the last line of the results file, after the 4th comma and before the 5th
Bstart=$(tail -n 1 "$resultsFile" | sed -n "s/^\([0-9.-]\+\),\([0-9.-]\+\),\([0-9.-]\+\),\([0-9.-]\+\),\([0-9.-]\+\),.*/\5/p");

# echoing Bstart
echo "Bstart = $Bstart"

# defining Bend as 10% higher
#Bend=$(bc <<< "scale=8;$Bstart*1.1");
Bend=$(python -c "print $Bstart*1.1");

# echoing B_end
echo "Bend = $Bend"

# writing B_start and B_end to inputs file
sed -i "s/^B\([^0-9.-]\+\)\([0-9.-]\+\) /B\1$Bstart/" $inputsFile
sed -i "s/^B\([^0-9.-]\+\)\([0-9.-]\+\)\([^0-9.-]\+\)\([0-9.-]\+\) /B\1\2\3$Bend/" $inputsFile

# submitting job
echo "submitting job"
msub $jobfile

echo "done - thank you lois xxx"
