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
file4="results/nr/nrmain_cosmos_5.dat"
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

#file4t u 6:($3==11? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, E=0.0" with points pointtype 1 lc rgb "red", \
	#file4t u 6:($3==11? ($5==1? ($8==0? $15: 1/0): 1/0): 1/0) title "K_g, E=0.0, " with points pointtype 1 lc rgb "orange", \

#set xrange [0.1:0.2]
#plot file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.02? $15: 1/0): 1/0): 1/0): 1/0) title "K_g, E=0.02, {\Symbol m}=1" with points #pointtype 1 lc rgb "orange", \
#	file4 u 6:($3==11? ($5==1? ($7==0.15? ($8==0.02? $15: 1/0): 1/0): 1/0): 1/0) title "K_g, E=0.04, {\Symbol m}=0.15" with points #pointtype 1 lc rgb "green", \
#	file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.04? $15: 1/0): 1/0): 1/0): 1/0) title "K_g, E=0.04, {\Symbol m}=1" with points #pointtype 1 lc rgb "violet", \
#	file4 u 6:($3==11? ($5==1? ($7==0.15? ($8==0.04? $15: 1/0): 1/0): 1/0): 1/0) title "K_g, E=0.04, {\Symbol m}=0.15" with points #pointtype 1 lc rgb "grey"

#set xrange [0.1:0.2]	
plot file4 u 6:($3==11? ($5==1? ($7==1? ($8==0? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=0" with points pointtype 1 lc rgb "violet" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.1? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=0.1, {\Symbol m}=1" with points pointtype 1 lc rgb "red" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.2? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=0.2, {\Symbol m}=1" with points pointtype 1 lc rgb "orange" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.3? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=0.3, {\Symbol m}=1" with points pointtype 1 lc rgb "blue" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.4? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=0.4, {\Symbol m}=1" with points pointtype 1 lc rgb "black" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==0.5? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=0.5, {\Symbol m}=1" with points pointtype 1 lc rgb "yellow" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==1? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=1.0, {\Symbol m}=1" with points pointtype 1 lc rgb "green" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==1.2? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=1.2, {\Symbol m}=1" with points pointtype 1 lc rgb "cyan" axes x1y1, \
	file5 u 6:($3==11? ($5==1? ($7==1? ($8==1.2? $9: 1/0): 1/0): 1/0): 1/0) notitle with points pointtype 1 lc rgb "cyan" axes x1y1, \
	file4 u 6:($3==11? ($5==1? ($7==1? ($8==1.4? $9: 1/0): 1/0): 1/0): 1/0) title "S, E=1.4, {\Symbol m}=1" with points pointtype 1 lc rgb "brown" axes x1y1, \
	file5 u 6:($3==11? ($5==1? ($7==1? ($8==1.4? $9: 1/0): 1/0): 1/0): 1/0) notitle with points pointtype 1 lc rgb "brown" axes x1y1, \


#set xrange [0.1:0.2]
#plot file4 u 6:($3==8? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, K=8" with points pointtype 1 lc rgb "green", \
#	file4 u 6:($3==8? ($5==1? ($8==0? $13: 1/0): 1/0): 1/0) title "K_g, K=8" with points pointtype 1 lc rgb "cyan", \
#	file4 u 6:($3==11? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, K=11" with points pointtype 1 lc rgb "red", \
#	file4 u 6:($3==11? ($5==1? ($8==0? $13: 1/0): 1/0): 1/0) title "K_g, K=11, " with points pointtype 1 lc rgb "orange", \
#	file4t u 6:($3==12? ($5==1? ($8==0? $9: 1/0): 1/0): 1/0) title "S, K=12" with points pointtype 1 lc rgb "blue", \
#	file4t u 6:($3==12? ($5==1? ($8==0? $13: 1/0): 1/0): 1/0) title "K_g, K=12" with points pointtype 1 lc rgb "violet", \
#	am(1.0,1.0) with lines title "Affleck-Manton result" lc rgb "black"
	
#plot file1 u ($3==11? ($5==1? ($6==0.15? $7: 1/0): 1/0): 1/0):($3==11? ($5==1? ($6==0.15? $10: 1/0): 1/0): 1/0) title "a=0.15, {/#Symbol m}=1.0" with points pointtype 1 lc rgb "red", \
#	file2 u ($3==11? ($5==1? ($6==0.15? $7: 1/0): 1/0): 1/0):($3==11? ($5==1? ($6==0.15? $8: 1/0): 1/0): 1/0) notitle with points #pointtype 1 lc rgb "red", \
#	file4 u ($3==11? ($6==0.15? ($7==0.15? $8: 1/0): 1/0): 1/0):($3==11? ($6==0.15? ($7==0.15? $9: 1/0): 1/0): 1/0) title "a=0.15, #{/Symbol m}=0.15" with points pointtype 2 lc rgb "blue", \
#	file4 u ($3==11? ($6==0.15? ($7==0.5? $8: 1/0): 1/0): 1/0):($3==11? ($6==0.15? ($7==0.5? $9: 1/0): 1/0): 1/0) title "a=0.15, {/#Symbol m}=0.5" with points pointtype 6 lc rgb "green", \
#	p(1.0,1.0,x) with lines title "weak coupling result" lc rgb "orange"

pause -1
