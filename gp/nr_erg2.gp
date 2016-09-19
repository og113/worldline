# second gnuplot program to plot erg conservation stuff


file="data/temp/160905171821P4_pl_0_run_2.dat"

E=0.000670932

set key below
set autoscale
set xtic auto
set ytic auto
set title "E vs s, fit with f=a+b*sin(2*{/Symbol p}*j/N)"
set xlabel "s"
set ylabel "E"

# fitting function
pi=3.1415926535897932;
K=11;
N=2.0**K;
w=(2.0*pi/N);
offset=0.01
f(x)= a + b*sin(w*(0.5+x)) + offset; # the offset helps the fitting process for some reason

gaussian=0
justTotal=0

# fitting data
fit f(x) file using 0:1 via a, b;
a1=a;
b1=b;

# titles
titleErg(t0,a0,b0) = sprintf("%s, a=%g, b=%g",t0,a0,b0)


plot file using 0:1 with lines lt 1 lc rgb "blue" title titleErg("Total",a1+offset,b1), \
	E with lines lt 1 lc rgb "red" title "dS/d{/ Symbol b}"

pause -1

