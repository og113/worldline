# gnuplot program to plot some early results of the temp calcs

file1="db/nr/plotsHand/s_len_dx_a_kappa_1_K_11_pot_9_T_0.02.dat"
file2="db/nr/plotsHand/s_len_dx_a_kappa_1_K_11_pot_8_T_0.02.dat"

set key below
set xtic auto
set ytic auto
set y2tic auto
#set logscale y
#set format y "10^{%L}"
set autoscale
set title "...; a=0.003, K=11, pot=12"
set xlabel "a"
set ylabel "S"
set y2label "dx/a"

plot file1 using 1:4 with points title "S, pot=9" lt 1 lc rgb "red", \
	file1 using 1:($2/$1) with points title "dx, pot=9" lt 1 lc rgb "blue", \
#	file2 using 1:4 with points title "S, pot=8" lt 1 lc rgb "green", \
#	file2 using 1:($2/$1) with points title "dx, pot=8" lt 1 lc rgb "orange", \

pause -1
