#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file1="data/temp/loop_G_1_B_1_M_1.435_a_0.15_mu_1.dat
file2="data/temp/loop_G_1_B_1_M_1.435_a_0.15_mu_0.15.dat
file3="data/temp/loop_G_1_B_0.98_M_1.439_a_0.15.dat"
file4="data/temp/loop_G_1_B_0.97_M_1.439_a_0.15.dat"
file5="data/temp/loop_G_1_B_0.96_M_1.439_a_0.15.dat"
file6="data/temp/loop_G_1_B_0.95_M_1.439_a_0.15.dat"
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
set title "Lemons, a=0.15, B=1.0, g=1.0, K=11, E=1.435"
set xlabel "z"
set ylabel "t"
set xrange [-0.02:0.02]
set yrange [1.3:1.4]
plot file1 using 3:4 title "{/Symbol m}=1.0" with lines lc rgb "red", \
	file2 using 3:4 title "{/Symbol m}=0.15" lt 2 lc rgb "blue"#, \
#	file3 using 3:4 title "B=0.98" with lines lc rgb "yellow", \
#	file4 using 3:4 title "B=0.97" with lines lc rgb "green", \
#	file5 using 3:4 title "B=0.96" with lines lc rgb "blue", \
#	file6 using 3:4 title "B=0.95" with lines lc rgb "violet"#, \
#	file7 using 3:4 title "B=0.94" with lines lc rgb "brown", \
#	file8 using 3:4 title "B=0.93" with lines lc rgb "black", \
#	file9 using 3:4 title "B=0.92" with lines lc rgb "grey", \
#	file10 using 3:4 title "B=0.91" with lines lc rgb "pink", \

pause -1

