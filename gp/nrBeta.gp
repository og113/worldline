#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# files
file2="results/nr/nr2.dat"
file4="results/nr/nr4.dat"

# approximate analytic results
pi=3.1415926535897932
AMscaled(kappa)= pi - kappa/4.0
wcscaled(E)= pi-2.0*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )
ptemp(T,N)=real(system(sprintf("./perturbativeFiniteTemp %f %f", T, N)))+0
phightemp(eta,kappa)= 2.0*(1.0-sqrt(kappa/4.0/pi))/eta
pweak(eta)= eta<0.5? pi: pi-2.0*( acos(1.0/2.0/eta) + (1.0/2.0/eta)*sqrt(1.0-(1.0/2.0/eta)**2) )

unset log
unset label
set key below
set xtic auto
set autoscale
#set logscale y
#set format y "10^{%L}"
set xlabel "T=1/{/Symbol b}"
set ylabel "S"

Napprox=100
K=10
pot=8
kappa=0.01

# before mu was a choice, mu=1
ST2l='file2 u 11:($3==pot? ($4==K? ($7=='
ST2m='? ($8=='
ST2r='? $12: 1/0): 1/0): 1/0): 1/0)'
ST4l='file4 u 11:($3==11? ($4==11? ($7=='
ST4r='? $12: 1/0): 1/0): 1/0)'

set xrange [0:2]
set yrange [1:4]
plot @ST2l 0.1 @ST2m 0.1 @ST2r t "Finite temperature calculation, closed topology, a=0.1, K=10" with points lc rgb "blue", \
	@ST4l 0.005 @ST4r t "Finite temperature calculation, open topology, a=0.005, K=11" with points lc rgb "magenta", \
	ptemp(x,Napprox) title "low temperature result" with lines lc rgb "red", \
	pweak(x) title "weak coupling result" with lines lc rgb "green", \
	phightemp(x,kappa) title "high temperature result" with lines lc rgb "orange"

pause -1
