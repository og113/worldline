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
unset key
set autoscale
set xtic auto
set ytic auto
set title "4d loop"

set multiplot layout 3,2 rowsfirst

set xlabel "x"
set ylabel "y"
plot file using 1:2 with lines

set xlabel "x"
set ylabel "z"
plot file using 1:3 with lines

set xlabel "x"
set ylabel "t"
plot file using 1:4 with lines

set xlabel "y"
set ylabel "z"
plot file using 2:3 with lines

set xlabel "y"
set ylabel "t"
plot file using 2:4 with lines

set xlabel "z"
set ylabel "t"
plot file using 3:4 with lines

unset multiplot

pause -1

