#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file="data/temp/xCurvature_K_11_G_1_B_1_M_1.43_run_1.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

stats file using 0 nooutput
N = STATS_max

stats file using 3:4 nooutput
min = (STATS_min_y < STATS_min_x ? STATS_min_y : STATS_min_x)
max = (STATS_max_y > STATS_max_x ? STATS_max_y : STATS_max_x)

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "4d loop"
set xlabel "z"
set ylabel "t"
#set xrange [min:max]
#set yrange [min:max]
#set size ratio -1
set xrange [-0.02:0.02]
set yrange [1.23:1.24]
#set size ratio -1

plot file using 3:4 with lines lc rgb 'black', \
	file using 3:4 ev N::0::N with lines lc rgb 'black' noti, \
	file using 3:4:6 with points pointtype 7 palette
	
pause -1

