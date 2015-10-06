#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
file="results/VvsK_dim_4.dat"

unset log
unset label
set autoscale
set xtic auto
set ytic auto
set title "v vs log_2(N)"
set xlabel "log_2(N)"
set ylabel "v"
f(x) = a+b*x
fit f(x) file using 3:11 via a, b
set key below
plot file using 3:11:12 with yerrorbars, \
	f(x) title sprintf("f(x) = %3.4f+%3.4f*x",a,b) with lines

pause -1
