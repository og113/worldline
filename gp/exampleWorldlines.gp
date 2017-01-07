#gnuplot program to plot some example worldline

#set termoption enhanced
#set encoding utf8 # for tilde

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

file1="data/temp/loop_kappa_0.4_E_0_a_0.02_T_0.5_pot_13.dat"
file2="data/temp/loop_kappa_0.4_E_0_a_0.02_T_0.62_pot_15.dat"
file3="data/temp/loop_kappa_0.4_E_0_a_0.02_T_0.7_pot_15.dat"
file4="data/temp/loop_kappa_0.5_E_0_a_0.02_T_0.71_pot_15.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

#set terminal latex
#set terminal latex
#set out 'temp/exampleInstantons.tex'

unset log
unset label
unset key
unset title
set autoscale
set xtic auto
set ytic auto
set xtics font "Arial,16" 
set ytics font "Arial,16" 
set xlabel font "Arial,20" 
set ylabel font "Arial,20" 


set multiplot layout 1,4 rowsfirst

set xlabel "x^3"
set ylabel "x^4" rotate by 0

set xrange [-1:1]
set yrange [-1:1]

set label 11 center at graph 0.5,char 1 "a)" font "Arial,24"
set bmargin 5
plot file1 using 1:2 with lines lt 1 lc rgb "black"

set label 11 "b)"
plot file2 using 1:2 with lines lt 1 lc rgb "black"

set label 11 "c)"
plot file3 using 1:2 with lines lt 1 lc rgb "black"

set label 11 "d)"
plot file4 using 1:2 with lines lt 1 lc rgb "black"

unset multiplot

#set out

pause -1

