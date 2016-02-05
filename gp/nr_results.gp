#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files where results are
file1="results/nr/nrmain_cosmos.dat"
file2="results/nr/nrmain_cosmos_2.dat"
file3="results/nr/nrmain_cosmos_3.dat"

# approximate analytic result
pi=3.1415926535897932
p(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

# fitting data
#f(x) = a+b*x
#fit f(x) file using (log($4)):(log($9)) via a, b

unset log
unset label
set key below
set xtic auto
#set logscale y
#set format y "10^{%L}"
set autoscale
set title "Induced instanton, g=1.0, B=1.0, K=10, E=1.4"
set xlabel "a"
set ylabel "S"

set xrange [0.13:0.15]
plot file3 u 6:($1==160203084722? $11: "1/0") title "NR result" ls 1 lc rgb "blue"#, \
	#file2 u 7:10 title "NR result, a=0.15" ls 1 lc rgb "green", \
	#file1 u 7:(p(1.0,1.0,$7)) title "weak coupling result" with lines lc rgb "red"

pause -1
