# v quick plot program for output of gsl_fft

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

unset log
unset label
set key below
set autoscale
set xtic auto
set ytic auto
set title "gsl fft test plot"
plot "data/temp/gsl_fft_x.dat" u 1:2 w lines lc rgb "red" t "data", \
	"data/temp/gsl_fft_xft.dat" u 1:2 w lines lc rgb "blue" t "fft of data"
	
pause -1
