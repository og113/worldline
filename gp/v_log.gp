#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \


unset log
unset label
set autoscale
set xtic auto
set ytic auto
set title "v vs log(N)"
set xlabel "log(N)"
set ylabel "v"
f(x) = a+b*x
fit f(x) "results/loop_dim_4.dat" using (log(2**$3)):8 via a, b
set key below
plot "results/loop_dim_4.dat" using (log(2**$3)):8 with lines, \
	f(x) title sprintf("f(x) = %3.4f+%3.4f*x",a,b) with lines

pause -1
