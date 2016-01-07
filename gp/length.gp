#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
file1="results/s0/length_s0_dim_4.dat"
file2="results/s0/length_s0+v_dim_4.dat"

unset log
unset label
#set logscale xy
set ytic auto
set xtic auto
set autoscale
#set format y "10^{%L}"
set title "log(length) vs log_2(N)"
set xlabel "log_2(N)"
set ylabel "log(length)"

# fitting data
f(x) = a+b*x
fit f(x) file1 using 3:(log($9)) via a, b
g(x) = c+d*x
fit g(x) file2 using 3:(log($9)) via c, d

set key below
set xrange [4:10]

plot file1 using 3:(log($9)):(log($10))with yerrorbars, \
	f(x) title sprintf("f(x) = %3.4f+%3.4f*x",a,b) with lines lc rgb "blue", \
	file2 using 3:(log($9)):(log($10)) with yerrorbars, \
	g(x) title sprintf("g(x) = %3.4f+%3.4f*x",c,d) with lines lc rgb "red"

pause -1
