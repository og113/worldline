# gnuplot program to plot z versus a

file1="db/nr/plots/Z_a_pot_11_K_10_kappa_0.01/z_a_pot_11_K_10_kappa_0.01.dat"
file2="data/nr/3dpot/extrema_bifurcate_0_pot_4_Kappa_0.01.dat"
file3="temp/zminmaxa.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Z_{max} and Z_{min} versus a"
set xlabel "a"
set ylabel "Z"

plot file1 using 1:2 with points ls 1  lt 4 lc rgb "orange" t "N-R z_{max}, from nr3", \
	file3 using 1:2 with points ls 1  lt 1 lc rgb "red" t "N-R z_{min}", \
	file3 using 1:3 with points ls 1  lt 1 lc rgb "black" t "N-R z_{max}", \
	file2 using 3:5 with points ls 1 lt 2 lc rgb "blue" t "3d minima", \
	file2 using 3:4 with points ls 1 lt 2 lc rgb "green" t "3d maxima", \

pause -1

