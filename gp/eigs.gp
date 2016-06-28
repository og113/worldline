# gnuplot program to plot eigenvalues

file="temp/eigs.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Eigenvalues"
set xlabel "n"
set ylabel "{/Symbol l}_n"

# fitting
pi=3.1415926535897932
beta=0.5
kappa=1
a0=2.0**8
a=((2.0*pi/beta)**2)/a0
b=2.0*kappa/16.0/pi
print a
print b
lambda(x)=a*(floor((x)/2)**2)-b
#fit lambda(x) file using 0:1 via a,b

plot file using 0:1 with points ls 1 noenhanced title "Numerical results", \
	lambda(x) with lines ls 1 lt 1 lc rgb "red" title sprintf("f(x) = %3.4f*x^2+%3.4f",a,b);

pause -1

