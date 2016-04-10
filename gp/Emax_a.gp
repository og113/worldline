#gnuplot program to plot Emax vs a from 3d potential calculation and from nr failure

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# approximate analytic result

# files
file="results/nr/Emax_Vmin_a.dat"

# threshold energy, from 3d analysis without cut off
pi=3.1415926535897932
E_threshold=2.0-1.0/sqrt(pi)

set key below
set autoscale
set title "E_{max} and V_{min} vs a, pot=4, {/Symbol l}=a"
set xlabel "a"
set ylabel "E"
set xrange [0:0.2]
set yrange [0:2]
set arrow from 0,E_threshold to 0.2,E_threshold nohead
set label "Threshold from 3d potential" at 0.1,(E_threshold+0.1)

plot file using 1:2:($2-0.00):($2+0.02) with yerrorbars title "E_{max} from 4d NR calculation" lc rgb "blue", \
	file u 1:3 title "V_{min} from 3d calculation" with points lc rgb "red"

pause -1
