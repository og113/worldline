#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \


unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "2d loop"

set multiplot layout 2,3 rowsfirst

set xlabel "x"
set ylabel "y"
plot "data/temp/projections.dat" using 1:2 with lines

set xlabel "x"
set ylabel "z"
plot "data/temp/projections.dat" using 2:3 with lines

set xlabel "x"
set ylabel "t"
plot "data/temp/projections.dat" using 1:3 with lines

set xlabel "y"
set ylabel "z"
plot "data/temp/projections.dat" using 3:4 with lines

set xlabel "y"
set ylabel "t"
plot "data/temp/projections.dat" using 3:4 with lines

set xlabel "z"
set ylabel "t"
plot "data/temp/projections.dat" using 3:4 with lines

pause -1

unset multiplot
