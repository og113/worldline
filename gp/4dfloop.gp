#
# This definition computes the sum of the first 10, 100, 1000 fourier
# coefficients of a (particular) square wave.

#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term fig size 1600,800; \
set output outFile; \

file="data/temp/floopAscii.dat"

fourier(a, k, t) = a*cos(2.0*pi*k*t)
sum(x,N) = sum [k=1:N] fourier(k, x)

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "floop"

set multiplot layout 2,2 rowsfirst

set xlabel "x"
set ylabel "y"
set zlabel "z"
splot file using 1:2:3 with lines

set xlabel "y"
set ylabel "z"
set zlabel "t"
splot file using 2:3:4 with lines

set xlabel "z"
set ylabel "t"
set zlabel "x"
splot file using 1:4 with lines

set xlabel "t"
set ylabel "x"
set zlabel "y"
splot file using 2:3 with lines

pause -1

unset multiplot












set title "Finite summation of 10, 100, 1000 fourier coefficients"

set samples 500
set xrange [-10:10]
set yrange [-0.4:1.2]
set key bottom right

fourier(k, x) = sin(3./2*k)/k * 2./3*cos(k*x)
sum10(x)   = 1./2 + sum [k=1:10]   fourier(k, x)
sum100(x)  = 1./2 + sum [k=1:100]  fourier(k, x)
sum1000(x) = 1./2 + sum [k=1:1000] fourier(k, x)

plot \
    sum10(x)   title "1./2 + sum [k=1:10]   sin(3./2*k)/k * 2./3*cos(k*x)", \
    sum100(x)  title "1./2 + sum [k=1:100]  sin(3./2*k)/k * 2./3*cos(k*x)", \
    sum1000(x) title "1./2 + sum [k=1:1000] sin(3./2*k)/k * 2./3*cos(k*x)"
