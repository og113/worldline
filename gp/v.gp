#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files for v0
#file1="data/s0/151016145509r_dim_4_K_10.dat"
#file2="data/s0/151016145518r_dim_4_K_10.dat"
#file3="data/s0/151016145527r_dim_4_K_10.dat"

# files for v1
file1="data/s0/151110162856v_dim_4_K_4.dat"
file2="data/s0/151110162900v_dim_4_K_6.dat"
file3="data/s0/151110162940v_dim_4_K_8.dat"

# finding mean
#mean1(x) = a1
#mean2(x) = a2
#mean3(x) = a3
#fit mean1(x) file1 using 9:10 via a1
#fit mean2(x) file2 using 9:10 via a2
#fit mean3(x) file3 using 9:10 via a3

#histogram stuff
n=50 #number of intervals
max=200.00 #max value
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
set title "v distribution, a=0.1"
set xlabel "v"
set ylabel "no. of groups"

#set xrange [0:200]    
plot file1 @hist1 title "K=4" lc rgb "red", \
	file2 @hist2 title "K=5" lc rgb "blue", \
	file3 @hist3 title "K=6" lc rgb "green"

pause -1
