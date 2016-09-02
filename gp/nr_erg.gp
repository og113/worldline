# gnuplot program to plot erg conservation stuff

file1="data/temp/160902121009P4_pl_0_run_4.dat"
file2="data/temp/160902121009P4_pl_1_run_3.dat"
file3="data/temp/160902121009P4_pl_2_run_3.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "E vs s, for different K"#fit with f=a+b*sin(2*{/Symbol p}*j/N)"
set xlabel "s"
set ylabel "E"

# fitting function
pi=3.1415926535897932;
K=10;
N=2.0**K;
w=(2.0*pi/N);
offset=0.01
f(x)= a + b*sin(w*(0.5+x)) + offset; # the offset helps the fitting process for some reason

gaussian=0
disjoint=0

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
fit f(x) file1 using 0:7 via a, b;
a7=a;
b7=b;

#if (gaussian==1) \
fit f(x) file1 using 0:6 via a, b;
a6=a;
b6=b;
if (gaussian==1) \
fit f(x) file1 using 0:8 via a, b;
a8=a;
b8=b;

if (gaussian==0) \
fit f(x) file1 using 0:5 via a, b;
a5=a;
b5=b;

sign(x,x0) = x<x0 ? 1 : -1;

# titles
titleErg(t0,a0,b0) = sprintf("%s, a=%g, b=%g",t0,a0,b0)

if (disjoint==1) \
plot file1 using 0:1 with lines lt 1 lc rgb "blue" title titleErg("Total",a1+offset,b1), \
	file1 using 0:2 with lines lt 1 lc rgb "red" title titleErg("Q S0",a2+offset,b2), \
	file1 using 0:3 with lines lt 1 lc rgb "green" title titleErg("Q I0",a3+offset,b3), \
	file1 using 0:4 with lines lt 1 lc rgb "orange" title titleErg("Q V",a4+offset,b4), \
	file1 using 0:6 with lines lt 1 lc rgb "brown" title titleErg("R Q",a7+offset,b7), \
	file1 using 0:7 with lines lt 1 lc rgb "cyan" title titleErg("Q G",a6+offset,b6), \
	file1 using 0:8 with lines lt 1 lc rgb "black" title titleErg("R G",a8+offset,b8), \
	#file1 using 0:5 with lines lt 1 lc rgb "magenta" title titleErg("Q L",a5+offset,b5), \
	
if (disjoint==0) \
plot file1 using ($0*4/256.0):1 with lines lt 1 lc rgb "blue" title "Total, K=6", \
	file2 using ($0*2/256.0):1 with lines lt 1 lc rgb "red" title "Total, K=7", \
	file3 using ($0/256.0):1 with lines lt 1 lc rgb "green" title "Total, K=8", \
	
	#titleErg("Total",a1+offset,b1), \
	
	#file1 using 0:2 with lines lt 1 lc rgb "red" title titleErg("Q S0",a2+offset,b2), \
	#file1 using 0:3 with lines lt 1 lc rgb "green" title titleErg("Q I0",a3+offset,b3), \
	#file1 using 0:4 with lines lt 1 lc rgb "orange" title titleErg("Q V",a4+offset,b4), \
	#file1 using 0:5 with lines lt 1 lc rgb "magenta" title titleErg("Q L",a5+offset,b5), \
	#file1 using 0:7 with lines lt 1 lc rgb "brown" title titleErg("R Q",a7+offset,b7), \
	#file1 using 0:6 with lines lt 1 lc rgb "cyan" title titleErg("Q G",a6+offset,b6), \
	#file1 using 0:8 with lines lt 1 lc rgb "black" title titleErg("R G",a8+offset,b8), \

# quick test of fit
#a7=a;
#b7=b;
#plot file1 using 0:7 with lines lt 1 lc rgb "blue" title titleErg("R Q",a7+offset,b7), \
#	f(x) with lines lt 1 lc rgb "red"

pause -1

