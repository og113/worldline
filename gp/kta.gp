# gnuplot program to plot parameter space explored for thermal calc

#file10="db/nr/plotsHand/kta_pot_10.dat"
#file11="db/nr/plotsHand/kta_pot_11.dat"
file12="db/nr/plotsHand/kta_pot_12.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Potential 11: Kappa, Temperature, a"
set xlabel "Kappa"
set ylabel "Temperature"
set zlabel "a"

splot file12 using 1:2:3 with points ls 1 lc rgb "blue" title "pot 12", \
#	file10 using 1:2:3 with points ls 1 lc rgb "green" title "pot 10", \
#	file11 using 1:2:3 with points ls 1 lc rgb "red" title "pot 11", \
	

pause -1

