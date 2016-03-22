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
filel7="results/nr/nrmain_laptop_7.dat"
filewkb="results/wkb/wkb_schwinger_rate.dat"

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
set ytic auto
set logscale y
set format y "10^{%L}"
set autoscale
set xlabel "E"
set ylabel "S"
set xtics 0.1
set grid xtics ytics

B=1
g=1
K=11

set title sprintf("Induced instanton, K=%g, {/Symbol m}=a",K)

ut6l='file6 u 9:($3==K? ($4==g? ($5==B? ($6==0? ($7=='
ut6m='? ($8=='
ut6r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0): 1/0)'
ut7l='file7 u 9:($3==K? ($4==g? ($5==B? ($6==0? ($7=='
ut7r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
utl7l='filel7 u 9:($3==K? ($4==g? ($5==B? ($6==0? ($7=='
utl7r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
style='with points'

plot @utl7l 0.15 @utl7r t "exp, a=0.15" @style pt 1 lc rgb "red", \
	@utl7l 0.125 @utl7r t "exp, a=0.125" @style pt 2 lc rgb "blue", \
	@utl7l 0.1 @utl7r t "exp, a=0.1" @style pt 4 lc rgb "green", \
	@ut7l 0.15 @ut7r t "original, a=0.15" @style pt 8 lc rgb "orange", \
	@ut7l 0.125 @ut7r t "original, a=0.125" @style pt 10 lc rgb "cyan", \
	@ut7l 0.1 @ut7r t "original, a=0.1" @style pt 1 lc rgb "pink", \
	@ut6l 0.15 @ut6m 0.15 @ut6r notitle @style pt 8 lc rgb "orange", \
	filewkb u 1:2 title "WKB non-relativistic result, a=0" @style pt 1 lc rgb "brown", \
	p(g,B,x) title "weak coupling result" lc rgb "black"


pause -1
