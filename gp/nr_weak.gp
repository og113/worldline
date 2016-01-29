#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# approximate analytic result
pi=3.1415926535897932
p(m,g,B)=pi/g/B - (2.0/g/B)*( asin(m/2.0) + (m/2.0)*sqrt(1-m*m/4.0) )

# files
file1="results/nr/nr_photon.dat"
file2="results/nr/nr_photon_nocuspterm.dat"

unset log
unset label
set key below
set logscale y
set xtic auto
set autoscale
set format y "10^{%L}"
set title "Induced NR results, K=8, g=1, B=1"
set xlabel "M"
set ylabel "S"

set xrange [0:2]
plot file1 u 7:10 title "NR result, a=0.15" ls 1 lc rgb "blue", \
	file2 u 7:8 title "NR result, no cusp regularisation term" ls 1 lc rgb "green", \
	file1 u 7:(p($7,1.0,1.0)) title "Weak coupling analytic result" with lines lc rgb "red"

pause -1
