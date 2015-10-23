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

# finding mean
mean1(x) = a1
mean2(x) = a2
mean3(x) = a3
fit mean1(x) file1 using 7:8 via a1
fit mean2(x) file2 using 7:8 via a2
fit mean3(x) file3 using 7:8 via a3

#histogram stuff
n=200 #number of intervals
max=10.00 #max value
min=0.00 #min value
width=(max-min)/n #interval width
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5
set macros
hist1='u (hist($8-a1,width)):(1.0/4096.0/width) smooth freq w boxes'
hist2='u (hist($8-a2,width)):(1.0/4096.0/width) smooth freq w boxes'
hist3='u (hist($8-a3,width)):(1.0/4096.0/width) smooth freq w boxes'

unset log
unset label
set key below
set autoscale
#set logscale y
set xtic auto
set ytic auto
set title "v distribution"
set xlabel "<v>"
set ylabel "no. of groups"

#set xrange [0:100]    
plot file1 @hist1 title "K=8" lc rgb "red", \
	file2 @hist2 title "K=10" lc rgb "blue", \
	file3 @hist3 title "K=12" lc rgb "green"

pause -1
