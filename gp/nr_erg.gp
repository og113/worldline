# gnuplot program to plot erg conservation stuff

file1="data/temp/160901132004P4_run_2.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "E vs j, fit with f=a+b*sin(2*{/Symbol p}*j/N)"
set xlabel "j"
set ylabel "E"

# fitting function
pi=3.1415926535897932;
K=10;
N=2.0**K;
w=(2.0*pi/N);
offset=0.01
f(x)= a + b*sin(w*(0.5+x)) + offset; # the offset helps the fitting process for some reason

# fitting data
fit f(x) file1 using 0:1 via a, b;
a1=a;
b1=b;
fit f(x) file1 using 0:2 via a, b;
a2=a;
b2=b;
fit f(x) file1 using 0:3 via a, b;
a3=a;
b3=b;
fit f(x) file1 using 0:4 via a, b;
a4=a;
b4=b;
#fit f(x) file1 using 0:5 via a, b;
#a5=a;
#b5=b;
fit f(x) file1 using 0:6 via a, b;
a6=a;
b6=b;
fit f(x) file1 using 0:7 via a, b;
a7=a;
b7=b;
fit f(x) file1 using 0:8 via a, b;
a8=a;
b8=b;

sign(x,x0) = x<x0 ? 1 : -1;

# titles
titleErg(t0,a0,b0) = sprintf("%s, a=%g, b=%g",t0,a0,b0)

disjoint=1

if (disjoint==1) \
plot file1 using 0:($1*sign($0,2**9)) with lines lt 1 lc rgb "blue" title titleErg("Total",a1+offset,b1), \
	file1 using 0:($2*sign($0,N/2)) with lines lt 1 lc rgb "red" title titleErg("Q S0",a2+offset,b2), \
	file1 using 0:($3*sign($0,N/2)) with lines lt 1 lc rgb "green" title titleErg("Q I0",a3+offset,b3), \
	file1 using 0:($4*sign($0,N/2)) with lines lt 1 lc rgb "orange" title titleErg("Q V",a4+offset,b4), \
	file1 using 0:($6*sign($0,N/2)) with lines lt 1 lc rgb "brown" title titleErg("R Q",a7+offset,b7), \
	file1 using 0:($7*sign($0,N/2)) with lines lt 1 lc rgb "cyan" title titleErg("Q G",a6+offset,b6), \
	file1 using 0:($8*sign($0,N/2)) with lines lt 1 lc rgb "black" title titleErg("R G",a8+offset,b8), \
	#file1 using 0:($5*sign($0,N/2)) with lines lt 1 lc rgb "magenta" title titleErg("Q L",a5+offset,b5), \
	
if (disjoint==0) \
plot file1 using 0:1 with lines lt 1 lc rgb "blue" title titleErg("Total",a1+offset,b1), \
	file1 using 0:2 with lines lt 1 lc rgb "red" title titleErg("Q S0",a2+offset,b2), \
	file1 using 0:3 with lines lt 1 lc rgb "green" title titleErg("Q I0",a3+offset,b3), \
	file1 using 0:4 with lines lt 1 lc rgb "orange" title titleErg("Q V",a4+offset,b4), \
	file1 using 0:5 with lines lt 1 lc rgb "magenta" title titleErg("Q L",a5+offset,b5), \
	file1 using 0:7 with lines lt 1 lc rgb "brown" title titleErg("R Q",a7+offset,b7), \
	#file1 using 0:6 with lines lt 1 lc rgb "cyan" title titleErg("Q G",a6+offset,b6), \
	#file1 using 0:8 with lines lt 1 lc rgb "black" title titleErg("R G",a8+offset,b8), \

# quick test of fit
#a7=a;
#b7=b;
#plot file1 using 0:7 with lines lt 1 lc rgb "blue" title titleErg("R Q",a7+offset,b7), \
#	f(x) with lines lt 1 lc rgb "red"

pause -1

