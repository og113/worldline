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
set title "Vacuum NR results, a=0.15, g=1.0, B=1.0, {\Symbol m}=1.0"
set xlabel "E"
set ylabel "S"
#set arrow from 0.5,6 to 0,(pi-0.25)
#set label "Affleck-Manton result" at 0.5,6
#set arrow from (2-1/sqrt(pi)),10 to (2-1/sqrt(pi)),0.1 nohead
#set label "2-sqrt(g^3*B/PI)" at (2-1/sqrt(pi)+0.05),1

set xrange [0:2]
plot file4 u 8:($3==11? ($5==1? ($6==0.15? ($7==1? $9: 1/0): 1/0): 1/0): 1/0) title "S, K=11" with points pointtype 1 lc rgb "blue", \
	file2 u 7:($3==12? ($5==1? ($6==0.15? ($4==1? $8: 1/0): 1/0): 1/0): 1/0) title "S, K=12" with points pointtype 2 lc rgb "green", \
	wc(1.0,1.0,x) title "Weak coupling result" with lines lc rgb "red"

pause -1

file2 u 5:7:($6==0.15? $8: "1/0") title "NR results, 2" ls 1 lc rgb "red", \
