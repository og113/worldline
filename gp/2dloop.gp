#gnuplot program to plot 2d worldline

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

file="data/temp/loopAscii.txt"

#if you want to choose which file to load
if (exists("inFile")) \
file=inFile; \

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "2d loop"
set xlabel "x"
set ylabel "y"
plot file using 1:2 with lines, \
	"data/temp/loopAsciiMet_run_1.dat" using 1:2 with lines, \
	"data/temp/loopAsciiMet_run_5.dat" using 1:2 with lines

pause -1
