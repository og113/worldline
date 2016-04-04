#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
fileb="data/nr/3dpot/extrema_bifurcate_1_Kappa_1.dat"
filenb="data/nr/3dpot/extrema_bifurcate_0_Kappa_1.dat"

unset log
unset label
set ytic auto
set xtic auto
set autoscale
set title "r vs a"
set xlabel "a"
set ylabel "r"
set key below

set xrange [0.01:0.1]

plot filenb using 1:2 title "r_{max}" with lines lc rgb "red", \
	filenb using 1:3 title "r_{min}" with lines lc rgb "blue"

pause -1
