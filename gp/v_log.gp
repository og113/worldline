#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
file="results/s0/v1VsK_a_0.dat"

unset log
unset label
set autoscale
set xtic auto
set ytic auto
set title "v1 vs log_2(N), a=0"
set xlabel "log_2(N)"
set ylabel "v"
f(x) = a+b*x
fit f(x) file using 3:14 via a, b
set key below
set xrange [2:7] 
plot file using 3:14:15 with yerrorbars, \
	f(x) title sprintf("f(x) = %3.4f+%3.4f*x",a,b) with lines

pause -1
