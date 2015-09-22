#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# approximate analytic result
p(x)=x/sinh(x)

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "<W>(TeB)"
set xlabel "TeB"
set ylabel "<W>

set xrange [0.0000001:5]    
plot "results/loop_dim_4.dat" u 6:9 ls 1, \
	p(x) with lines

pause -1
