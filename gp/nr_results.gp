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
set logscale y
set format y "10^{%L}"
set autoscale
set title "Induced instanton, g=1.0, B=1.0, K=11"
set xlabel "E"
set ylabel "S"

set xrange [0:1.5]
plot file1 u ($3==11? ($5==1? ($6==0.15? $7: 1/0): 1/0): 1/0):($3==11? ($5==1? ($6==0.15? $10: 1/0): 1/0): 1/0) title "a=0.15, {/Symbol m}=1.0" with points pointtype 1 lc rgb "red", \
	file2 u ($3==11? ($5==1? ($6==0.15? $7: 1/0): 1/0): 1/0):($3==11? ($5==1? ($6==0.15? $8: 1/0): 1/0): 1/0) notitle with points pointtype 1 lc rgb "red", \
	file4 u ($3==11? ($6==0.15? ($7==0.15? $8: 1/0): 1/0): 1/0):($3==11? ($6==0.15? ($7==0.15? $9: 1/0): 1/0): 1/0) title "a=0.15, {/Symbol m}=0.15" with points pointtype 2 lc rgb "blue", \
	p(1.0,1.0,x) with lines title "weak coupling result" lc rgb "black"

pause -1
