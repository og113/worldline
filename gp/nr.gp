#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# files
file1="results/nr/nrmain_cosmos_a_0.1.dat"
file24="results/nr/nrmain_cosmos_24format.dat"
file2="results/nr/nrmain_cosmos_2.dat"
file3="results/nr/nrmain_cosmos_3.dat"
file4="results/nr/nrmain_cosmos_4.dat"
file5="results/nr/nrmain_cosmos_5.dat"
file6="results/nr/nrmain_cosmos_6.dat"

# approximate analytic results
pi=3.1415926535897932
AM(g,B)=pi/g/B - g**2/4.0
wc(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

unset log
unset label
set key below
set logscale y
set xtic auto
set autoscale
set format y "10^{%L}"
set xlabel "E"
set ylabel "S"
#set arrow from 0.5,6 to 0,(pi-0.25)
#set label "Affleck-Manton result" at 0.5,6
#set arrow from (2-1/sqrt(pi)),10 to (2-1/sqrt(pi)),0.001 nohead
#set label "2-sqrt(g^3*B/PI)" at (2-1/sqrt(pi)+0.05),1

a=0.15
B=1
g=1
K=12
Ng=0

SE2='file2 u 7:($3==K? ($5==B? ($6==a? $8: 1/0): 1/0): 1/0)'
SE3='file3 u 7:($3==K? ($5==B? ($6==a?$8: 1/0): 1/0): 1/0)'
SE4l='file4 u 8:($3==K? ($5==B? ($6==a? ($7=='
SE4r='? $9: 1/0): 1/0): 1/0): 1/0)'
SE5l='file5 u 9:($3==K? ($5==B? ($6==Ng? ($7==a? ($8=='
SE5r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
SE6l='file6 u 9:($3==K? ($5==B? ($6==Ng? ($7==a? ($8=='
SE6r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
style1='with points pointtype 1 lc rgb "red"'
style2='with points pointtype 2 lc rgb "blue"'
titleS(mu) = sprintf("{/Symbol m}=%g",mu)
set title sprintf("Vacuum NR results, a=%g, g=%g, B=%g, K=%g",a,g,B,K)

set xrange [0:2]
plot @SE4l 0.15 @SE4r t titleS(0.15) @style1, \
	@SE5l 0.15 @SE5r notitle @style1, \
	@SE6l 0.15 @SE6r notitle @style1, \
	@SE2 t titleS(1) @style2, \
	@SE3 notitle @style2, \
	@SE4l 1 @SE4r notitle @style2, \
	@SE5l 1 @SE5r notitle @style2, \
	@SE6l 1 @SE6r notitle @style2, \
	#wc(1.0,1.0,x) title "Weak coupling result" with lines lc rgb "red"

pause -1

file2 u 5:7:($6==0.15? $8: "1/0") title "NR results, 2" ls 1 lc rgb "red", \
