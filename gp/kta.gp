# gnuplot program to plot parameter space explored for thermal calc

#file10="db/nr/plotsHand/kta_pot_10.dat"
#file11="db/nr/plotsHand/kta_pot_11.dat"
file12="db/nr/plotsHand/kta_pot_12.dat"
file8="db/nr/plotsHand/kta_pot_8.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Potential 12: Kappa, Temperature, a"
set xlabel "{/Symbol k}"
set ylabel "Temperature"
set zlabel "a"

pi=3.1415926535897932
Tmax(kappa)=sqrt(2)/pi**0.75/kappa**0.25

splot file12 using 1:2:3 with points ls 1 lc rgb "blue" title "pot 12", \
	file12 using 1:(Tmax($1)):(0.01) with lines ls 1 lc rgb "black" title "T_{Max}", \
#	file10 using 1:2:3 with points ls 1 lc rgb "green" title "pot 10", \
#	file11 using 1:2:3 with points ls 1 lc rgb "red" title "pot 11", \
	

pause -1

