#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files where results are
file0="results/nr/nrmain.dat"
file1="results/nr/nrmain_cosmos.dat"
file2="results/nr/nrmain_cosmos_2.dat"
file3="results/nr/nrmain_cosmos_3.dat"
file4="results/nr/nrmain_cosmos_4.dat"
file5="results/nr/nrmain_cosmos_5.dat"
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
set title "Instantons, g=1.0, B=1.0, K=K"
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
K=11
mu=1

ut4l='file4 u 6:($3==K? ($5==B? ($7==mu? ($8=='
ut4r='? $9: 1/0): 1/0): 1/0): 1/0)'
ut5l='file4 u 7:($3==K? ($5==B? ($8==mu? ($9=='
ut5r='? $10: 1/0): 1/0): 1/0): 1/0)'
style='with points pointtype 1'
titleS(E) = sprintf("S, E=%g, {\Symbol m}=%g",E,mu)

set xrange [0.01:1]	
plot @ut4l 0 @ut4r t titleS(0) @style lc rgb "red", \
	@ut5l 0 @ut5r notitle @style lc rgb "red", \
	@ut4l 0.2 @ut4r t titleS(0.2) @style lc rgb "orange", \
	@ut5l 0.2 @ut5r notitle @style lc rgb "orange", \
	@ut4l 0.3 @ut4r t titleS(0.3) @style lc rgb "green", \
	@ut5l 0.3 @ut5r notitle @style lc rgb "green", \
	@ut4l 0.4 @ut4r t titleS(0.4) @style lc rgb "blue", \
	@ut5l 0.4 @ut5r notitle @style lc rgb "blue", \
	@ut4l 0.5 @ut4r t titleS(0.5) @style lc rgb "black", \
	@ut5l 0.5 @ut5r notitle @style lc rgb "black", \
	@ut4l 1.0 @ut4r t titleS(1.0) @style lc rgb "yellow", \
	@ut5l 1.0 @ut5r notitle @style lc rgb "yellow", \
	@ut4l 1.2 @ut4r t titleS(1.2) @style lc rgb "cyan", \
	@ut5l 1.2 @ut5r notitle @style lc rgb "cyan", \
	@ut4l 1.4 @ut4r t titleS(1.4) @style lc rgb "brown", \
	@ut5l 1.4 @ut5r notitle @style lc rgb "orange", \


#set xrange [0.1:0.2]
#plot file4 u 6:($3==8? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, K=8" with points pointtype 1 lc rgb "green", \
#	file4 u 6:($3==8? ($5==1? ($8==0? $13: 1/0): 1/0): 1/0) title "K_g, K=8" with points pointtype 1 lc rgb "cyan", \
#	file4 u 6:($3==K? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, K=K" with points pointtype 1 lc rgb "red", \
#	file4 u 6:($3==K? ($5==1? ($8==0? $13: 1/0): 1/0): 1/0) title "K_g, K=K, " with points pointtype 1 lc rgb "orange", \
#	file4t u 6:($3==12? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, K=12" with points pointtype 1 lc rgb "blue", \
#	file4t u 6:($3==12? ($5==1? ($8==0? $13: 1/0): 1/0): 1/0) title "K_g, K=12" with points pointtype 1 lc rgb "violet", \
#	am(1.0,1.0) with lines title "Affleck-Manton result" lc rgb "black"
	
#plot file1 u ($3==K? ($5==1? ($6==0.15? $7: 1/0): 1/0): 1/0):($3==K? ($5==1? ($6==0.15? $10: 1/0): 1/0): 1/0) title "a=0.15, {/#Symbol m}=1.0" with points pointtype 1 lc rgb "red", \
#	file2 u ($3==K? ($5==1? ($6==0.15? $7: 1/0): 1/0): 1/0):($3==K? ($5==1? ($6==0.15? $8: 1/0): 1/0): 1/0) notitle with points #pointtype 1 lc rgb "red", \
#	file4 u ($3==K? ($6==0.15? ($7==0.15? $8: 1/0): 1/0): 1/0):($3==K? ($6==0.15? ($7==0.15? $9: 1/0): 1/0): 1/0) title "a=0.15, #{/Symbol m}=0.15" with points pointtype 2 lc rgb "blue", \
#	file4 u ($3==K? ($6==0.15? ($7==0.5? $8: 1/0): 1/0): 1/0):($3==K? ($6==0.15? ($7==0.5? $9: 1/0): 1/0): 1/0) title "a=0.15, {/#Symbol m}=0.5" with points pointtype 6 lc rgb "green", \
#	p(1.0,1.0,x) with lines title "weak coupling result" lc rgb "orange"

pause -1
