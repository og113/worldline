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
set ztic auto
set xtics 0.5
set ytics 0.5
set ztics 0.5
set title "4d loop"

set multiplot layout 2,2 rowsfirst

set xlabel "x"
set ylabel "y"
set zlabel "z"
splot file using 1:2:3 with lines

set xlabel "y"
set ylabel "z"
set zlabel "t"
splot file using 2:3:4 with lines

set xlabel "z"
set ylabel "t"
set zlabel "x"
splot file using 3:4:1 with lines

set xlabel "t"
set ylabel "x"
set zlabel "y"
splot file using 4:1:2 with lines

pause -1

unset multiplot
