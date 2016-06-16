#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term svg size 1600,800; \
set output outFile; \

# files
file="results/nr/nr_error3.dat"

unset log
unset label
set key below
set xtic auto
set ytic auto
set autoscale
set xlabel "T"
set ylabel "a"#{/Symbol k}"
set zlabel "Sol"
set title "Errors of failed N-R calculations"

plot file using 11:7:13 with points pointtype 7 palette

pause -1
