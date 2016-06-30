#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \


file="data/temp/160630114523mdsEarly1_K_9_kappa_0.01_E_0_a_0.01_mu_0.1_pot_11_T_1.8_run_1.dat"
file2="data/temp/160630134229mdsEarly1_K_9_kappa_0.01_E_0_a_0.01_mu_0.1_pot_17_T_1.8_run_1.dat"
file3="data/temp/160630103011mdsEarly1_K_9_kappa_0.01_E_0_a_0.01_mu_0.1_pot_17_T_1.8_run_1.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

offset=0.0002
offset4=0.0002

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "4d loop"
set xlabel "point"
set ylabel "loop coordinate"
plot file using 0:($3+offset) t "z" with points pt 1 lc rgb "blue", \
	file using 0:($4+offset4) t "t" with points pt 1 lc rgb "green", \
	file using 0:($1+offset) t "x" with points pt 1 lc rgb "red", \
	file using 0:($2+offset) t "y" with points pt 1 lc rgb "black", \
#	file2 using 0:($3+offset) t "z2" with points pt 1 lc rgb "orange", \
#	file2 using 0:($4+offset4) t "t2" with points pt 1 lc rgb "cyan", \
#	file2 using 0:($1+offset) t "x2" with points pt 1 lc rgb "magenta", \
#	file2 using 0:($2+offset) t "y2" with points pt 1 lc rgb "brown", \

pause -1

