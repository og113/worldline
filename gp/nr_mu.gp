#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files where results are
file0="results/nr/nrmain.dat"
file1="results/nr/nrmain_cosmos.dat"
file2="results/nr/nrmain_cosmos_2.dat"
file3="results/nr/nrmain_cosmos_3.dat"
file4="results/nr/nrmain_cosmos_4.dat"
file5="results/nr/nrmain_cosmos_5.dat"
file6="results/nr/nrmain_cosmos_6.dat"
file7="results/nr/nrmain_cosmos_7.dat"
file4t="temp/nrmain_cosmos_4.dat"
file5t="temp/nrmain_cosmos_5.dat"

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
set logscale y
set format y "10^{%L}"
set autoscale
set xlabel "{/Symbol l}"
set ylabel "angle"
set grid xtics ytics

a=0.15
B=1
E=0.02
mu=0.15

set title sprintf("Induced instanton, a=%g, E=%g",a,E)

ut7l='file7 u 8:($9==E? ($5==B? ($6==0? ($7==a? ($3=='
ut7r='? -$19: 1/0): 1/0): 1/0): 1/0): 1/0)'
style='with points'
titleS(K) = sprintf("angle neighbouring cusp, K=%g",K)

set xrange [0.1:0.2]	
plot @ut7l 8 @ut7r t titleS(8) @style pt 2 lc rgb "black", \
	@ut7l 9 @ut7r t titleS(9) @style pt 2 lc rgb "orange", \
	@ut7l 10 @ut7r t titleS(10) @style pt 2 lc rgb "red", \
	@ut7l 11 @ut7r t titleS(11) @style pt 4 lc rgb "green", \
	@ut7l 12 @ut7r t titleS(12) @style pt 8 lc rgb "blue", \
	#@ut7l 1.2 @ut7r t titleS(1.2) @style pt 10 lc rgb "violet", \
	

pause -1
