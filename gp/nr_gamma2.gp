#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
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
file4t="temp/nrmain_cosmos_4.dat"
file5t="temp/nrmain_cosmos_5.dat"

# approximate analytic result
pi=3.1415926535897932
#p(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

# fitting data
#f(x) = a+b*x
#fit f(x) file using (log($4)):(log($9)) via a, b

# gamma free
mass(g,a)=(1.0-g**2/8.0/pi/a)

# approximate analytic results
pi=3.1415926535897932
mass(g,a)=1-g**2/8.0/pi/a
gamma_free(E,g,a)=2.0*asin(E/2.0/mass(g,a))

#set logscale y
#set format y "10^{%L}"
set xtic auto
set autoscale
set xlabel "E"
set ylabel "{/Symbol g}"

B=1
g=1
K=11
Ng=0

# 
GE4l='file4 u 8:($3==K? ($5==B? ($6=='
GE4m='? ($7=='
GE4r='? $10: 1/0): 1/0): 1/0): 1/0)'
GE5l='file5 u 9:($3==K? ($5==B? ($6==Ng? ($7=='
GE5m='? ($8=='
GE5r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
GE6l='file6 u 9:($3==K? ($5==B? ($6==Ng? ($7=='
GE6m='? ($8=='
GE6r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
GE7l='file7 u 9:($3==K? ($5==B? ($6==Ng? ($7=='
GE7m='? ($8=='
GE7r='? $11: 1/0): 1/0): 1/0): 1/0): 1/0)'
style1='with points pointtype 1 lc rgb "red"'
style2='with points pointtype 2 lc rgb "blue"'
style3='with points pointtype 4 lc rgb "green"'
style4='with points pointtype 8 lc rgb "violet"'
style5='with points pointtype 10 lc rgb "orange"'
style6='with points pointtype 19 lc rgb "cyan"'
titleS(a,mu) = sprintf("a=%g, {/Symbol m}=%g",a,mu)
titleF(a)=sprintf("Free result, a=%g",a)
set title sprintf("NR results, g=%g, B=%g, K=%g",g,B,K)

set xrange [0:2]
plot @GE4l 0.2 @GE4m 0.2 @GE4r t titleS(0.2,0.2) @style4, \
	@GE5l 0.2 @GE5m 0.2 @GE5r notitle @style4, \
	@GE6l 0.2 @GE6m 0.2 @GE6r notitle @style4, \
	@GE7l 0.2 @GE7m 0.2 @GE7r notitle @style4, \
	gamma_free(x,g,0.2) title titleF(0.2) with lines lc rgb "violet", \
	@GE4l 0.15 @GE4m 0.15 @GE4r t titleS(0.15,0.15) @style1, \
	@GE5l 0.15 @GE5m 0.15 @GE5r notitle @style1, \
	@GE6l 0.15 @GE6m 0.15 @GE6r notitle @style1, \
	@GE7l 0.15 @GE7m 0.15 @GE7r notitle @style1, \
	gamma_free(x,g,0.15) title titleF(0.15) with lines lc rgb "red", \
	@GE4l 0.125 @GE4m 0.125 @GE4r t titleS(0.125,0.125) @style3, \
	@GE5l 0.125 @GE5m 0.125 @GE5r notitle @style3, \
	@GE6l 0.125 @GE6m 0.125 @GE6r notitle @style3, \
	@GE7l 0.125 @GE7m 0.125 @GE7r notitle @style3, \
	gamma_free(x,g,0.125) title titleF(0.125) with lines lc rgb "green", \
	@GE4l 0.1 @GE4m 0.1 @GE4r t titleS(0.1,0.1) @style2, \
	@GE5l 0.1 @GE5m 0.1 @GE5r notitle @style2, \
	@GE6l 0.1 @GE6m 0.1 @GE6r notitle @style2, \
	@GE7l 0.1 @GE7m 0.1 @GE7r notitle @style2, \
	gamma_free(x,g,0.1) title titleF(0.1) with lines lc rgb "blue", \
	@GE4l 0.075 @GE4m 0.075 @GE4r t titleS(0.075,0.075) @style5, \
	@GE5l 0.075 @GE5m 0.075 @GE5r notitle @style5, \
	@GE6l 0.075 @GE6m 0.075 @GE6r notitle @style5, \
	@GE7l 0.075 @GE7m 0.075 @GE7r notitle @style5, \
	gamma_free(x,g,0.075) title titleF(0.075) with lines lc rgb "orange", \
	@GE4l 0.05 @GE4m 0.05 @GE4r t titleS(0.05,0.05) @style6, \
	@GE5l 0.05 @GE5m 0.05 @GE5r notitle @style6, \
	@GE6l 0.05 @GE6m 0.05 @GE6r notitle @style6, \
	@GE7l 0.05 @GE7m 0.05 @GE7r notitle @style6, \
	gamma_free(x,g,0.05) title titleF(0.05) with lines lc rgb "cyan", \
	gamma_free(x,0,0.15) title "Free result without counterterms" with lines lc rgb "black"
	
pause -1
