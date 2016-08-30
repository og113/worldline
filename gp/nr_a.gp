# gnuplot program to plot a dependence

file1="db/nr/plotsHand/s_a_dx_T_0.02_pot_9_K_11_kappa_1.dat"
file2="db/nr/plotsHand/s_a_dx_T_0.2_pot_9_K_11_kappa_1.dat"
file3="db/nr/plotsHand/s_a_dx_T_0.4_pot_9_K_11_kappa_1.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "S vs a; {/Symbol k}=1, K=11, pot=9"
set xlabel "a"
set ylabel "S"
set cblabel "dx/a"

set yrange [2.8:2.92]

plot file3 using 1:2:($3/$1) with points lt 1 palette title "S, dx/a, T=0.4", \
	file2 using 1:2:($3/$1) with points lt 2 palette title "S, dx/a, T=0.2", \
	file1 using 1:2:($3/$1) with points lt 4 palette title "S, dx/a, T=0.02"

pause -1

