# gnuplot program to plot eigenvalues

file="data/temp/160621163708eigenvalues_pl_0_run_1.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Eigenvalues"
set xlabel "n"
set ylabel "{/Symbol l}_n"

# fitting
pi=3.1415926535897932
beta=0.8
kappa=0.0193765
#a=2.0*((2.0*pi/beta)**2)
#b=4.0*sqrt(pi/kappa)*beta
lambda(x)=a*(floor(x/2+1)**2)+b
fit lambda(x) file using 0:1 via a,b

print lambda(0)

plot file using 0:1 with lines ls 1 noenhanced title "Numerical results", \
	lambda(x) with lines ls 1 lt 1 lc rgb "red" title sprintf("f(x) = %3.4f*x^2+%3.4f",a,b);

pause -1

