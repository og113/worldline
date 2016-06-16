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
set ylabel "{/Symbol k}"
set zlabel "a"
set title "Errors of failed N-R calculations"

splot file using 11:5:7:15 title "dSdx norm" with points pointtype 7 palette

pause -1
