#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
file="data/temp/s0.dat"

unset log
set key below
set autoscale
set xtic auto
set ytic auto
#set logscale y
#set format y "10^{%L}"
set title "V1 metropolis run, S=S0+V1, pi*L/a regularisation, Epsi=0.1, initially circle R=1"
set xlabel "t"
set ylabel "V1"
set key below
plot file using 0:1 with lines lc rgb "blue" title "correlator(t)"

pause -1
