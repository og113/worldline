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
g=1.0
B=1.0
maximum(g,B)=sqrt(g/4.0/pi/B) # 1.771395468764214?
minimum=0.189673005548958 # for a=0.05

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

if (exists("max")) \
set arrow from (-maximum(g,B)/2.0),0 to (maximum(g,B)/2.0),0 nohead

if (exists("min")) \
set arrow from (-minimum/2.0),1 to (minimum/2.0),1 nohead

plot file using 3:4 with lines

pause -1

