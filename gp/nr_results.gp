#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("oSEFile")) \
if (oSEFile ne 'gui') \
set term png size 1600,800; \
set oSEpSE oSEFile; \

# files where results are
file0="results/nr/nrmain.dat"
file1="results/nr/nrmain_cosmos.dat"
file2="results/nr/nrmain_cosmos_2.dat"
file3="results/nr/nrmain_cosmos_3.dat"
file4="results/nr/nrmain_cosmos_4.dat"
file5="results/nr/nrmain_cosmos_5.dat"
file6="results/nr/nrmain_cosmos_6.dat"
file4t="temp/nrmain_cosmos_4.dat"
file5t="temp/nrmain_cosmos_5.dat"

# approximate analytic result
pi=3.1415926535897932
am(g,B) = pi/g/B - g**2/4.0
p(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

# fitting data
#f(x) = a+b*x
#fit f(x) file using (log($4)):(log($9)) via a, b

unset log
unset label
set key below
set xtic auto
set logscale xy
set format y "10^{%L}"
set autoscale
set title "Instantons, g=1.0, B=1.0, K=11"
set xlabel "a"
set ylabel "S"

#set y2tic auto
#set logscale y2
#set y2label "K_g"
#set format y2 "10^{%L}"
#set y2label "K_g maximum"

# function for derivatives
d(x,y) = ($0 == 0) ? (x1 = x, y1 = y, 1/0) : (x2 = x1, x1 = x, y2 = y1, y1 = y, (y1-y2)/(x1-x2))

#set arrow from (2-1/sqrt(pi)),10 to (2-1/sqrt(pi)),0.1 nohead

a=0.15
B=1
Ng=0
K=11
mu=1

SE4l='file4 u 6:($3==K? ($5==B? ($7==mu? ($8=='
SE4r='? $9: 1/0): 1/0): 1/0): 1/0)'
SE5l='file5 u 7:($3==K? ($5==B? ($6==Ng? ($8==mu? ($9=='
SE5r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
SE6l='file6 u 7:($3==K? ($5==B? ($6==Ng? ($8==mu? ($9=='
SE6r='? $10: 1/0): 1/0): 1/0): 1/0): 1/0)'
style='with points pointtype 1'
titleS(E) = sprintf("S, E=%g, {/Symbol m}=%g",E,mu)

set xrange [0.01:1]	
plot @SE4l 0 @SE4r t titleS(0) @style lc rgb "red", \
	@SE5l 0 @SE5r notitle @style lc rgb "red", \
	@SE6l 0 @SE6r notitle @style lc rgb "red", \
	@SE4l 0.2 @SE4r t titleS(0.2) @style lc rgb "orange", \
	@SE5l 0.2 @SE5r notitle @style lc rgb "orange", \
	@SE6l 0.2 @SE6r notitle @style lc rgb "orange", \
	@SE4l 0.3 @SE4r t titleS(0.3) @style lc rgb "green", \
	@SE5l 0.3 @SE5r notitle @style lc rgb "green", \
	@SE6l 0.3 @SE6r notitle @style lc rgb "green", \
	@SE4l 0.4 @SE4r t titleS(0.4) @style lc rgb "blue", \
	@SE5l 0.4 @SE5r notitle @style lc rgb "blue", \
	@SE6l 0.4 @SE6r notitle @style lc rgb "blue", \
	@SE4l 0.5 @SE4r t titleS(0.5) @style lc rgb "black", \
	@SE5l 0.5 @SE5r notitle @style lc rgb "black", \
	@SE6l 0.5 @SE6r notitle @style lc rgb "black", \
	@SE4l 0.7 @SE4r t titleS(0.7) @style lc rgb "pink", \
	@SE5l 0.7 @SE5r notitle @style lc rgb "pink", \
	@SE6l 0.7 @SE6r notitle @style lc rgb "pink", \
	@SE4l 1.0 @SE4r t titleS(1.0) @style lc rgb "yellow", \
	@SE5l 1.0 @SE5r notitle @style lc rgb "yellow", \
	@SE6l 1.0 @SE6r notitle @style lc rgb "yellow", \
	@SE4l 1.2 @SE4r t titleS(1.2) @style lc rgb "cyan", \
	@SE5l 1.2 @SE5r notitle @style lc rgb "cyan", \
	@SE6l 1.2 @SE6r notitle @style lc rgb "cyan", \
	@SE4l 1.4 @SE4r t titleS(1.4) @style lc rgb "brown", \
	@SE5l 1.4 @SE5r notitle @style lc rgb "brown", \
	@SE6l 1.4 @SE6r notitle @style lc rgb "brown", \

pause -1
