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
file8="results/nr/nrmain_cosmos_8.dat"

# approximate analytic result
pi=3.1415926535897932
p(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

# perturbative low temp result


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
set xrange[0:2]

B=1
g=1
K=11
pot=6

set title sprintf("Original regularisation with Gaussian repulsion, K=%g, {/Symbol m}=a",K)

ut6l='file6 u 9:($3==K? ($4==g? ($5==B? ($6==0? ($7=='
ut6m='? ($8=='
ut6r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0): 1/0)'
ut7l='file7 u 9:($3==K? ($4==g? ($5==B? ($6==0? ($7=='
ut7r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
utl7l='filel7 u 9:($3==K? ($4==g? ($5==B? ($6==0? ($7=='
utl7r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
ut8l='file8 u 11:($3=='
ut8m='? ($4==K? ($5==g? ($6==B? ($7==0? ($8=='
ut8r='? $12: 1/0): 1/0): 1/0): 1/0): 1/0): 1/0)'
style='with points'

if (pot==0) \
plot @ut7l 0.15 @ut7r t "original, a=0.15" @style pt 1 lc rgb "red", \
	@ut6l 0.15 @ut6m 0.15 @ut6r notitle @style pt 1 lc rgb "red", \
	@ut7l 0.125 @ut7r t "original, a=0.125" @style pt 2 lc rgb "blue", \
	@ut7l 0.1 @ut7r t "original, a=0.1" @style pt 4 lc rgb "green", \
	filewkb u 1:2 title "WKB non-relativistic result, a=0" with lines lc rgb "cyan", \
	p(g,B,x) title "weak coupling result" lc rgb "black"

if (pot==2) \
plot @utl7l 0.15 @utl7r t "exp, a=0.15" @style pt 1 lc rgb "red", \
	@utl7l 0.125 @utl7r t "exp, a=0.125" @style pt 2 lc rgb "blue", \
	@utl7l 0.1 @utl7r t "exp, a=0.1" @style pt 4 lc rgb "green", \
	filewkb u 1:2 title "WKB non-relativistic result, a=0" with lines lc rgb "cyan", \
	p(g,B,x) title "weak coupling result" lc rgb "black"
	
if (pot==4) \
plot @ut8l 4 @ut8m 0.15 @ut8r t "original + gaussian, a=0.15" @style pt 6 lc rgb "orange", \
	@ut8l 4 @ut8m 0.125 @ut8r t "original + gaussian, a=0.125" @style pt 8 lc rgb "cyan", \
	@ut8l 4 @ut8m 0.1 @ut8r t "original + gaussian, a=0.1" @style pt 10 lc rgb "pink", \
	@ut8l 4 @ut8m 0.075 @ut8r t "original + gaussian, a=0.075" @style pt 1 lc rgb "blue", \
	@ut8l 4 @ut8m 0.05 @ut8r t "original + gaussian, a=0.05" @style pt 2 lc rgb "red", \
	filewkb u 1:2 title "WKB non-relativistic result, a=0" with lines lc rgb "green", \
	p(g,B,x) title "weak coupling result" lc rgb "black"

if (pot==6) \
plot @ut8l 6 @ut8m 0.15 @ut8r t "exp + gaussian, a=0.15" @style pt 6 lc rgb "orange", \
	@ut8l 6 @ut8m 0.125 @ut8r t "exp + gaussian, a=0.125" @style pt 8 lc rgb "cyan", \
	@ut8l 6 @ut8m 0.1 @ut8r t "exp + gaussian, a=0.1" @style pt 10 lc rgb "pink", \
	@ut8l 6 @ut8m 0.075 @ut8r t "exp + gaussian, a=0.075" @style pt 1 lc rgb "blue", \
	@ut8l 6 @ut8m 0.05 @ut8r t "exp + gaussian, a=0.05" @style pt 2 lc rgb "red", \
	filewkb u 1:2 title "WKB non-relativistic result, a=0" with lines lc rgb "green", \
	p(g,B,x) title "weak coupling result" lc rgb "black"

pause -1
