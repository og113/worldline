#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
file="data/temp/corr.dat"

unset log
set key below
set autoscale
set xtic auto
set ytic auto
#set logscale y
#set format y "10^{%L}"
set title "correlator(t)"
set xlabel "t"
set ylabel "correlator"
set key below
plot file using 0:1 with lines lc rgb "blue" title "correlator(t)"

pause -1
