# gnuplot program to plot a dependence

file1="db/nr/plotsHand/s_a_dx_K_12_pot_9_T_0.02.dat"
file2="db/nr/plotsHand/s_a_dx_K_11_pot_9_T_0.02.dat"

file3="db/nr/plotsHand/s_a_dx_K_12_pot_9_T_0.4.dat"
file4="db/nr/plotsHand/s_a_dx_K_11_pot_9_T_0.4.dat"

file5="db/nr/plotsHand/s_a_dx_K_12_pot_8_T_0.02.dat"
file6="db/nr/plotsHand/s_a_dx_K_11_pot_8_T_0.02.dat"
file7="db/nr/plotsHand/s_a_dx_K_10_pot_8_T_0.02.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "S vs a; {/Symbol k}=1, pot=8"
set xlabel "a"
set ylabel "S"
set cblabel "dx/a"

y(x)=a*x+b
fit [0.02:0.1] y(x) file5 using 1:3 via a,b

#set yrange [2.8:2.86]

plot file5 using 1:3:($2/$1) with points lt 1 palette title "S, dx/a, K=12, T=0.02", \
	file6 using 1:3:($2/$1) with points lt 2 palette title "S, dx/a, K=11, T=0.02", \
	file7 using 1:3:($2/$1) with points lt 4 palette title "S, dx/a, K=11, T=0.02", \
	y(x) with lines ls 1 lt 1 lc rgb "red" title sprintf("y(x) = %3.4f*x+%3.4f",a,b);

pause -1

