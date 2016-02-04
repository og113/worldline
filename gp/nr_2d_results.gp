#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file where results are
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
set logscale z
set format z "10^{%L}"
set autoscale
set title "Induced instanton, g=1.0, a=0.15, K=11"
set xlabel "B"
set ylabel "E"
set zlabel "S"

set xrange [0:1]
set yrange [0:1.8]
splot file1 u 5:7:($6==0.15? ($3==11?$10:"1/0"): "1/0") title "NR results, 1" ls 1 lc rgb "blue", \
	file2 u 5:7:($6==0.15? $8: "1/0") title "NR results, 2" ls 1 lc rgb "red", \
	file3 u 5:7:($6==0.15? $8: "1/0") title "NR results, 3" ls 1 lc rgb "violet", \
	p(1.0,x,y) title "weak coupling result" with lines lc rgb "green"

pause -1
