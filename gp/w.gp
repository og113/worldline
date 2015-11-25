#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# approximate analytic result
p(x)=(x/sinh(x)-1-x**2/6.0)

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "<W>(TeB)"
set xlabel "TeB"
set ylabel "<W>

set xrange [0:10]    
plot "results/s0/floopOffice_dim_4.dat" u 8:(($11-1-$8**2/6.0)):12 with yerrorbars ls 1, \
	p(x) with lines

pause -1
