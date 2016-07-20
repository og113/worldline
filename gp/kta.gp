# gnuplot program to plot parameter space explored for thermal calc

file10="db/nr/plotsHand/pot_10_kta.dat"
file11="db/nr/plotsHand/pot_11_kta.dat"
file12="db/nr/plotsHand/pot_12_kta.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Potential 11: Kappa, Temperature, a"
set xlabel "Kappa"
set ylabel "Temperature"
set zlabel "a"

splot file10 using 1:2:3 with points ls 1 lc rgb "green" title "pot 10", \
	file11 using 1:2:3 with points ls 1 lc rgb "red" title "pot 11", \
	file12 using 1:2:3 with points ls 1 lc rgb "blue" title "pot 12", \

pause -1

