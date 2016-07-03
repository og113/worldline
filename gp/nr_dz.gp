# gnuplot program to plot dsdz_RmL and delta_RmL

file1="results/nr/nr_dz_pot_17.dat"
file2="results/nr/nr_dz_pot_15.dat"
file3="results/nr/nr_dz_pot_11.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set y2tic auto
set title "'Force' versus 'dz'; {/Symbol k}=0.01, a=0.005"
set xlabel "z"
set ylabel "dVdz, the 'Force'"
set y2label "dz"


plot file1 using 1:2 with lines lt 1 lc rgb "red" title "nonrel potential: dVdz", \
	file1 using 1:3 with lines lt 1 lc rgb "orange" title "nonrel potential: dz", \
	file2 using 1:2 with lines lt 1 lc rgb "blue" title "external potential: dVdz", \
	file2 using 1:3 with lines lt 1 lc rgb "black" title "external potential: dz", \
	file3 using 1:2 with lines lt 1 lc rgb "green" title "full thermal potential: dVdz", \
	file3 using 1:3 with lines lt 1 lc rgb "brown" title "full thermal potential: dz", \
	0 with lines	

pause -1

