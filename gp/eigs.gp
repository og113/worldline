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
beta=0.5555555555555
kappa=0.01
a0=2.0**8
a=((2.0*pi/beta)**2)/4.0/a0
b=4.0*sqrt(pi/kappa)*beta/a0
#print a
#print b
lambda(x)=a*(floor((x+1)/2)**2)-b
fit lambda(x) file using 0:1 via b

plot file using 0:1 with points ls 1 noenhanced title "Numerical results", \
	lambda(x) with lines ls 1 lt 1 lc rgb "red" title sprintf("f(x) = %3.4f*x^2+%3.4f",a,b);

pause -1

