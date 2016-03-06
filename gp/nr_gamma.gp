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

# approximate analytic results
pi=3.1415926535897932
mass(g,a)=1-g**2/8.0/pi/a
gamma_free(E,g,a)=2.0*asin(E/2.0/mass(g,a))

unset log
unset label
set key below
#set logscale y
#set format y "10^{%L}"
set xtic auto
set autoscale
set xlabel "E"
set ylabel "{/Symbol g}"

B=1
g=1
K=12
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
style2='with points pointtype 8 lc rgb "blue"'
titleS(a,mu) = sprintf("a=%g, {/Symbol m}=%g",a,mu)
set title sprintf("NR results, g=%g, B=%g, K=%g",g,B,K)

set xrange [0:2]
plot @GE4l 0.15 @GE4m 0.15 @GE4r t titleS(0.15,0.15) @style1, \
	@GE5l 0.15 @GE5m 0.15 @GE5r notitle @style1, \
	@GE6l 0.15 @GE6m 0.15 @GE6r notitle @style1, \
	@GE7l 0.15 @GE7m 0.15 @GE7r notitle @style1, \
	@GE4l 0.1 @GE4m 0.1 @GE4r t titleS(0.1,0.1) @style2, \
	@GE5l 0.1 @GE5m 0.1 @GE5r notitle @style2, \
	@GE6l 0.1 @GE6m 0.1 @GE6r notitle @style2, \
	@GE7l 0.1 @GE7m 0.1 @GE7r notitle @style2, \
	gamma_free(x,g,0.15) title "Free result including counterterms, a=0.15" with lines lc rgb "red", \
	gamma_free(x,g,0.1) title "Free result including counterterms, a=0.1" with lines lc rgb "blue", \
	gamma_free(x,0,0.15) title "Free result without counterterms" with lines lc rgb "black"
	

pause -1
