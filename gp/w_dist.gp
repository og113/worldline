#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# file
file1="data/s0/151016145509r_dim_4_K_10.dat"
file2="data/s0/151016145518r_dim_4_K_10.dat"
file3="data/s0/151016145527r_dim_4_K_10.dat"

#histogram stuff
n=40 #number of intervals
max=1.00 #max value
min=0.00 #min value
width=(max-min)/n #interval width
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5
set macros
hist1='u (hist($10,width)):(1.0/64.0/width) smooth freq w boxes'
hist2='u (hist($10,width)):(1.0/64.0/width) smooth freq w boxes'
hist3='u (hist($10,width)):(1.0/64.0/width) smooth freq w boxes'

unset log
unset label
set key below
set autoscale
#set logscale y
set xtic auto
set ytic auto
set title "w distribution"
set xlabel "<w>"
set ylabel "no. of groups"

#set xrange [0:100]    
plot file1 @hist1 title "B=0" lc rgb "red", \
	file2 @hist2 title "B=1" lc rgb "blue", \
	file3 @hist3 title "B=2" lc rgb "green"

pause -1
