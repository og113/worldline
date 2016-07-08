# gnuplot program to plot some early results of the disjoint temp calcs

file1="temp/firstResults_T_a_0.07_kappa_0.933055_K_11.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "z_{max} vs T;  {/Symbol k}=0.933055, a=0.07, K=11"
set xlabel "T"
set ylabel "z_{max}"


plot file1 using 1:3 with linespoints lt 1 lc rgb "red"

pause -1
