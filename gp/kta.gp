# gnuplot program to plot parameter space explored for thermal calc

#file10="db/nr/plotsHand/kta_pot_10.dat"
#file11="db/nr/plotsHand/kta_pot_11.dat"
file12="db/nr/plotsHand/kta_pot_12.dat"
file8="db/nr/plotsHand/kta_pot_8.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Potential 8: Kappa, Temperature, a"
set xlabel "{/Symbol k}"
set ylabel "Temperature"
set zlabel "a"

splot file8 using 1:2:3 with points ls 1 lc rgb "blue" title "pot 8", \
#	file10 using 1:2:3 with points ls 1 lc rgb "green" title "pot 10", \
#	file11 using 1:2:3 with points ls 1 lc rgb "red" title "pot 11", \
	

pause -1

