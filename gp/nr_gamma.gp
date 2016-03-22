#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files
file1="results/nr/nrmain_cosmos.dat"
file2="results/nr/nrmain_cosmos_2.dat"
file3="results/nr/nrmain_cosmos_3.dat"
file4="results/nr/nrmain_cosmos_4.dat"
file5="results/nr/nrmain_cosmos_5.dat"
file6="results/nr/nrmain_cosmos_6.dat"
file7="results/nr/nrmain_cosmos_7.dat"
filel7="results/nr/nrmain_laptop_7.dat"

# approximate analytic results
pi=3.1415926535897932
mass(g,a)=1-g**2/8.0/pi/a
gamma_free(E,g,a)=2.0*asin(E/2.0/mass(g,a))

unset log
unset label
set key below
set xtic auto
set ytic auto
#set logscale y
#set format y "10^{%L}"
set autoscale
set xlabel "E"
set ylabel "{/Symbol g}"
#set xtics 0.01
set grid xtics ytics

a=0.10
B=1
g=1
mu=0.10

set title sprintf("Induced instanton, a=%g, mu=%g",a,mu)

ut5l='file5 u 9:($7==a? ($5==B? ($6==0? ($8==mu? ($3=='
ut5r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
ut6l='file6 u 9:($8==mu? ($5==B? ($6==0? ($7==a? ($3=='
ut6r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
ut7l='file7 u 9:($8==mu? ($5==B? ($6==0? ($7==a? ($3=='
ut7r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
utl7l='filel7 u 9:($8==mu? ($5==B? ($6==0? ($7==a? ($3=='
utl7r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
style='with points'
titleS(K) = sprintf("{/Symbol g}, K=%g",K)

#set xrange [0.1:0.2]	
plot @utl7l 8 @utl7r t titleS(8) @style pt 1 lc rgb "red", \
	@utl7l 9 @utl7r t titleS(9) @style pt 2 lc rgb "orange", \
	@utl7l 10 @utl7r t titleS(10) @style pt 4 lc rgb "violet", \
	@utl7l 11 @utl7r t titleS(11) @style pt 4 lc rgb "green", \
	@utl7l 12 @utl7r t titleS(12) @style pt 8 lc rgb "blue", \
	@utl7l 13 @utl7r t titleS(13) @style pt 10 lc rgb "black", \
	@ut7l 8 @ut7r notitle @style pt 1 lc rgb "red", \
	@ut7l 9 @ut7r notitle @style pt 2 lc rgb "orange", \
	@ut7l 10 @ut7r notitle @style pt 4 lc rgb "violet", \
	@ut7l 11 @ut7r notitle @style pt 4 lc rgb "green", \
	@ut7l 12 @ut7r notitle @style pt 8 lc rgb "blue", \
	@ut7l 13 @ut7r notitle @style pt 10 lc rgb "black", \
	@ut6l 8 @ut6r notitle @style pt 1 lc rgb "red", \
	@ut6l 9 @ut6r notitle @style pt 2 lc rgb "orange", \
	@ut6l 10 @ut6r notitle @style pt 4 lc rgb "violet", \
	@ut6l 11 @ut6r notitle @style pt 4 lc rgb "green", \
	@ut6l 12 @ut6r notitle @style pt 8 lc rgb "blue", \
	@ut6l 13 @ut6r notitle @style pt 10 lc rgb "black", \
	@ut5l 8 @ut5r notitle @style pt 1 lc rgb "red", \
	@ut5l 9 @ut5r notitle @style pt 2 lc rgb "orange", \
	@ut5l 10 @ut5r notitle @style pt 4 lc rgb "violet", \
	@ut5l 11 @ut5r notitle @style pt 4 lc rgb "green", \
	@ut5l 12 @ut5r notitle @style pt 8 lc rgb "blue", \
	@ut5l 13 @ut5r notitle @style pt 10 lc rgb "black", \
	gamma_free(x,0.0,a) t "{/Symbol g}, continuum free theory" with lines lc rgb "black"

pause -1
