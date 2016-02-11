#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files where results are
file1="results/nr/nrmain_cosmos_4.dat"
file2="results/nr/nrmain_failures.dat"

# approximate analytic result
#pi=3.1415926535897932
#p(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

# fitting data
#f(x) = a+b*x
#fit f(x) file using (log($4)):(log($9)) via a, b

unset log
unset label
set key below
set xtic auto
set ytic auto
#set y2tic auto
#set logscale z
#set format z "10^{%L}"
set autoscale
set title "Induced instanton, a=0.15, K=11"
set xlabel "{/Symbol m}"
set ylabel "d(dx^2)/2dx^2"
set zlabel "E maximum"
#set y2label "average K*a"

set xrange [0:1.2]
#set yrange [0:1.8]
splot file2 u 7:15:8 title "E maximum" ls 1 lc rgb "blue"

#plot file1 u 7:($1==160208165923?$15:"1/0") title "d(dx^2)/2dx^2 at maximum E" ls 1 lc rgb "blue" axes x1y1, \
#	file1 u 7:($1==160208165923?$16:"1/0") title "average K*a" ls 1 lc rgb "red" axes x1y2

pause -1
