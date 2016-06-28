#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file="data/temp/151201162640mdsEarly1_K_8_G_1_B_1_run_1.dat"

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
set xlabel "point"
set ylabel "loop coordinate"
plot file using 0:3 t "z" with points pt 1 lc rgb "blue", \
	file using 0:4 t "t" with points pt 1 lc rgb "green", \
	file using 0:1 t "x" with points pt 1 lc rgb "red", \
	file using 0:2 t "y" with points pt 1 lc rgb "black", \

pause -1

