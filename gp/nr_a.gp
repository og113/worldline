# gnuplot program to plot a dependence

file1="db/nr/plotsHand/s_a_dx_K_12_pot_9_T_0.02.dat"
file2="db/nr/plotsHand/s_a_dx_K_11_pot_9_T_0.02.dat"

file3="db/nr/plotsHand/s_a_dx_K_12_pot_9_T_0.4.dat"
file4="db/nr/plotsHand/s_a_dx_K_11_pot_9_T_0.4.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "S vs a; {/Symbol k}=1, pot=9"
set xlabel "a"
set ylabel "S"
set cblabel "dx/a"

set yrange [2.8:2.86]

plot file3 using 1:3:($2/$1) with points lt 1 palette title "S, dx/a, K=12, T=0.4", \
	file4 using 1:3:($2/$1) with points lt 2 palette title "S, dx/a, K=11, T=0.4", \

pause -1

