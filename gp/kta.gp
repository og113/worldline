# gnuplot program to plot parameter space explored for thermal calc

file="results/nr/kta.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Kappa, Temperature, a"
set xlabel "Kappa"
set ylabel "Temperature"
set zlabel "a"

splot file using 1:2:3 with points ls 1 title "Parameters with curved solution"

pause -1

