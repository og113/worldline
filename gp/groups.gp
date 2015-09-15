#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

#histogram stuff
n=60 #number of intervals
max=60.00 #max value
min=0.00 #min value
width=(max-min)/n #interval width
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5
set macros
hist='u (hist($5,width)):(1.0/10000.0/width) smooth freq w boxes'

# approximate analytic result
pi=3.141592653589793
dof=31*4
c=-1+0.5*dof
p(x)=exp(-(x-c)+c*log(x/c))/sqrt(2.0*c*pi)

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "grouped results"
set xlabel "<S0>_g"
set ylabel "no. of groups"

set xrange [0:100]    
plot "results/150914115841loopGroups_dim_4_K_5.dat" @hist ls 1, \
	p(x) with lines

pause -1
