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
set title "Induced instanton, a=0.15, K=11, {/Symbol m}=0.05"
set xlabel "E"
set ylabel "cc max"

a=0.15
B=1
K=11
mu=0.15

ut6l='file6 u 7:($3==K? ($5==B? ($6==0? ($8==mu? ($9=='
ut6r='? $18: 1/0): 1/0): 1/0): 1/0): 1/0)'
style='with points'
titleS(E) = sprintf("cc max, E=%g, {/Symbol m}=%g",E,mu)

#set xrange [0.01:1]	
plot @ut6l 0.4 @ut6r t titleS(0.4) @style pt 2 lc rgb "orange", \
	@ut6l 0.2 @ut6r t titleS(0.2) @style pt 4 lc rgb "green", \
	@ut6l 0.7 @ut6r t titleS(0.7) @style pt 8 lc rgb "blue", \
	@ut6l 1.2 @ut6r t titleS(1.2) @style pt 10 lc rgb "violet", \
	

pause -1
