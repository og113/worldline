#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file1="data/temp/loop_kappa_0.1_E_0_a_0.01_T_0.02_pot_13.dat"
file2="data/temp/loop_kappa_0.1_E_0_a_0.01_T_0.02_pot_13.dat"
file3="data/temp/loop_kappa_0.1_E_0_a_0.01_T_0.02_pot_13.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "Finite temperature worldlines"

set multiplot layout 1,3 rowsfirst

set xlabel "x^3"
set ylabel "x^4"

plot file1 using 1:2 with lines lt 1 lc rgb "black" title "Worldline 1"
plot file2 using 1:2 with lines lt 1 lc rgb "black" title "Worldline 2"
plot file3 using 1:2 with lines lt 1 lc rgb "black" title "Worldline 3"

unset multiplot

pause -1

