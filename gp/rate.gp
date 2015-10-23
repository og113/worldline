#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# approximate analytic result
pi=3.1415926535897932
p(x)=real(system(sprintf("./schwingerRate %f", x)))+0

unset log
unset label
set key below
set logscale xy
set xtic auto
set autoscale
set format y "10^{%L}"
set title "rate(eE)"
set xlabel "eE"
set ylabel "rate"

set xrange [0.1:10]
plot "results/s0/Fr_dim_4.dat" u 7:($7**2*$11/16.0/pi**2):($12) with yerrorbars title "Monte-Carlo rate"ls 1 lc rgb "blue", \
	-p(x) title "Schwinger result" with lines lc rgb "red"

pause -1
