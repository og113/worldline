# gnuplot program to plot number of eigenvalues

file="results/nr/nr_eigs_pot_11_K_11.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set ztic auto
set title "Eigenvalues: pot=11; K=11; input=nonrel; run=1"
set xlabel "kappa"
set ylabel "T"
set zlabel "# eigenvalues"

splot file using 5:10:13 with points ls 1 lt 1 lc rgb "blue" title "negative eigenvalues", \
	file using 5:10:14 with points ls 1 lt 1 lc rgb "red" title "zero eigenvalues"

pause -1

