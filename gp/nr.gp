#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# approximate analytic result
pi=3.1415926535897932
p(g,B)=pi/g/B - g**2/4.0

# files
file1="results/nr/nr_vacuum_K.dat"
file2="results/nr/nr_vacuum_K.dat"

unset log
unset label
set key below
set logscale y
set xtic auto
set autoscale
set format y "10^{%L}"
set title "Vacuum NR results, a=0.1, g=1, B=1"
set xlabel "K"
set ylabel "S"

set xrange [5:11]
plot file u 3:9 title "NR result" ls 1 lc rgb "blue", \
	p(1.0,1.0) title "Analytic result" with lines lc rgb "red"

pause -1
