#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file0="data/temp/loop_G_1_B_1_M_0.04_a_0.13_mu_0.15.dat"
file1="data/temp/loop_G_1_B_1_M_0.04_a_0.14_mu_0.15.dat"
file2="data/temp/loop_G_1_B_1_M_0.04_a_0.15_mu_0.15.dat"
file3="data/temp/loop_G_1_B_1_M_0.04_a_0.16_mu_0.15.dat"
file4="data/temp/loop_G_1_B_1_M_0.04_a_0.17_mu_0.15.dat"
file5="data/temp/loop_G_1_B_1_M_1.42_a_0.15_mu_1_Ng_11.dat"
file6="data/temp/loop_G_1_B_1_M_1.42_a_0.15_mu_1_Ng_13.dat"
file7="data/temp/loop_G_1_B_0.4_M_1.439_a_0.15.dat"
file8="data/temp/loop_G_1_B_0.3_M_1.439_a_0.15.dat"
file9="data/temp/loop_G_1_B_0.2_M_1.439_a_0.15.dat"
file10="data/temp/loop_G_1_B_0.1_M_1.439_a_0.15.dat"


#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "Lemons, {/Symbol m}=0.15, B=1.0, g=1.0, K=11, E=0.04"
set xlabel "z"
set ylabel "t"
set xrange [-0.01:0.01]
set yrange [0.994:0.998]
style='with linespoints'
plot file0 using 3:4 title "a=0.13" @style lc rgb "black", \
	file1 using 3:4 title "a=0.14" @style lc rgb "red", \
	file2 using 3:4 title "a=0.15" @style lc rgb "blue", \
	file3 using 3:4 title "a=0.16" @style lc rgb "yellow", \
	file4 using 3:4 title "a=0.17" @style lc rgb "green", \
#	file5 using 3:4 title "Nc=11" @style lc rgb "violet", \
#	file6 using 3:4 title "B=0.95" @style lc rgb "violet"#, \
#	file7 using 3:4 title "B=0.94" @style lc rgb "brown", \
#	file8 using 3:4 title "B=0.93" @style lc rgb "black", \
#	file9 using 3:4 title "B=0.92" @style lc rgb "grey", \
#	file10 using 3:4 title "B=0.91" @style lc rgb "pink", \

pause -1

