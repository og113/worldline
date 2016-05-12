#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# files
file2="results/nr/nr2.dat"

# approximate analytic results
pi=3.1415926535897932
AM(g,B)=pi/g/B - g**2/4.0
wc(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )
ptemp(g,B,T)=AM(g,B)-g**2*((pi**4/90.0)*T**4+4.0*(pi**6/945.0)*T**6)
ptemp(T,N)=real(system(sprintf("./perturbativeFiniteTemp %f %f", T, N)))+0

unset log
unset label
set key below
set xtic auto
set autoscale
set logscale y
set format y "10^{%L}"
set xlabel "T"
set ylabel "S"

Napprox=100
K=10
pot=8

# before mu was a choice, mu=1
ST2l='file2 u 11:($3==pot? ($4==K? ($7=='
ST2m='? ($8=='
ST2r='? $12: 1/0): 1/0): 1/0): 1/0)'

set xrange [0:1.0]
plot @ST2l 0.1 @ST2m 0.1 @ST2r t "Finite temperature calculation" with points lc rgb "blue", \
	ptemp(x,Napprox) title "Perturbative result" with lines lc rgb "red"

pause -1
