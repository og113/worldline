#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file="data/temp/151201162640mdsEarly1_K_8_G_1_B_1_run_1.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "4d loop"
set xlabel "z"
set ylabel "t"
<<<<<<< HEAD
<<<<<<< HEAD
#set xrange [-0.04:0.04]
#set yrange [-0.99:-1.00]
plot file using 3:4 with linespoints pt 2
=======
set xrange [-1:1]
set yrange [-1:1]
=======
#set xrange [-1:1]
#set yrange [-1:1]
>>>>>>> 84ea179b219fcb515e65c74bbee74c0c93d487de
plot file using 3:4 with lines
>>>>>>> 13b684a57852489990d96a2f984e1aaf644e824f

pause -1

