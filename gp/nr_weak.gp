#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# approximate analytic result
pi=3.1415926535897932
p(m,g,B)=pi/g/B - (2.0/g/B)*( asin(m/2.0) + (m/2.0)*sqrt(1-m*m/4.0) )

# files
file1="results/nr/nr_a_0.15_K_11.dat"
#file2="results/nr/nr_a_0.10_K_10.dat"
#file3="results/nr/nr_a_0.05_K_10.dat"
#file4="results/nr/nr_a_0.02_K_10.dat"

unset log
unset label
set key below
set logscale y
set xtic auto
set autoscale
set format y "10^{%L}"
set title "Induced NR results, K=11, g=1, B=1, a=0.15"
set xlabel "E"
set ylabel "S"
#set arrow from 0.5,6 to 0,(pi-0.25)
#set label "Affleck-Manton result" at 0.5,6
#set arrow from (2-1/sqrt(pi)),10 to (2-1/sqrt(pi)),0.1 nohead
#set label "2-sqrt(g^3*B/PI)" at (2-1/sqrt(pi)+0.05),1

set xrange [0:2]
plot file1 u 7:10 title "NR result" ls 1 lc rgb "blue", \
	file1 u 7:(p($7,1.0,1.0)) title "Weak coupling analytic result" with lines lc rgb "red"
	#file2 u 7:10 title "NR result, a=0.10" ls 1 lc rgb "blue", \
	#file3 u 7:10 title "NR result, a=0.05" ls 1 lc rgb "green", \
	#file4 u 7:10 title "NR result, a=0.02" ls 1 lc rgb "red", \

pause -1
