#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file="data/temp/loop_G_1_B_1_M_1.22_a_0.05_mu_0.05_pot_4.dat"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

# stationary point of 3d potential
pi=3.1415926535897932
g=1
B=1
stationaryPoint(g,B)=sqrt(1/16.0/pi/g/B)

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "4d loop"
set xlabel "z"
set ylabel "t"
set xrange [-1:1]
#set yrange [-1:1]

if (exists("arrow")) \
set arrow from (-stationaryPoint(g,B)/2.0),0 to (stationaryPoint(g,B)/2.0),0 nohead


plot file using 3:4 with lines

pause -1

