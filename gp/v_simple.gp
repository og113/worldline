#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file in
file="data/s0/151109174721v_dim_4_K_10.dat"

#histogram stuff
n=10 #number of intervals
max=50.00 #max value
min=0.00 #min value
width=(max-min)/n #interval width
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5
set macros
hist='u (hist(-$10,width)):(1.0/64.0/width) smooth freq w boxes'

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "V1 distribution, K=6"
set xlabel "V1"
set ylabel "no. of groups"

set xrange [0:600]    
plot file @hist ls 1

pause -1
