#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \


file="data/temp/highTemp_kappa_0.310023_T_0.8_rank_0.dat"
file2="data/temp/addVectors.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "4d loop"
set xlabel "z"
set ylabel "t"
plot file using 3:4 t "original worldline" with points pt 1 lc rgb "blue", \
	file2 using 3:4 t "plus change" with points pt 1 lc rgb "green", \
	
pause -1
