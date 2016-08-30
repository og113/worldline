# gnuplot program to plot dsdz_RmL and delta_RmL

file1="db/nr/plotsHand/zmin_T_a_0.003_kappa_0.01_K_11_pot_12.dat"


set key below
set autoscale
set xtic auto
set ytic auto
set title "'dz' versus T; {/Symbol k}=0.01, a=0.003, K=11, pot=12"
set xlabel "Temp"
set ylabel "dz"

a=0.003

plot file1 using 1:2 with points lt 1 lc rgb "red" title "minimum dz", \
	a with lines lt 1 lc rgb "blue" title "a", \
	file1 using 1:($3*a) with points lt 1 lc rgb "black" title "dx"

pause -1

