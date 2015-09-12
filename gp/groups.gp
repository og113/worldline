#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

#histogram stuff
n=10 #number of intervals
max=20.00 #max value
min=0.00 #min value
width=(max-min)/n #interval width
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5
set macros
hist='u (hist($6,width)):(1.0) smooth freq w boxes'

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "grouped results"
set xlabel "<S0>_g"
set ylabel "no. of groups"
     
plot "results/150912085150loopGroups_dim_4_K_10.dat" @hist ls 1

pause -1
