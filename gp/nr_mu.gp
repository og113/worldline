#gnuplot program to plot results from each group of loops, should give a gaussian plot

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

# files where results are
file1="results/nr/nrmain_cosmos_4.dat"
file2="results/nr/nrmain_failures.dat"

# approximate analytic result
#pi=3.1415926535897932
#p(g,B,E)=pi/g/B-(2.0/g/B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-(E/2.0)**2) )

# fitting data
#f(x) = a+b*x
#fit f(x) file using (log($4)):(log($9)) via a, b

unset log
unset label
set key below
set xtic auto
set ytic auto
#set y2tic auto
#set logscale z
#set format z "10^{%L}"
set autoscale
set title "Induced instanton, a=0.15, K=11"
set xlabel "mu"
set ylabel "K_G*a maximum"
#set y2label "S"
#set zlabel "E maximum"
#set y2label "average K*a"

#set xrange [0:1.2]
#set yrange [0:1.8]
#splot file2 u 7:15:8 title "E maximum" ls 1 lc rgb "blue"

plot file1 u 7:($1==160208165923?($6==0.15?($8==1.43?$15:"1/0"):"1/0"):"1/0") title "NR results" ls 1 lc rgb "blue" axes x1y1

#plot file1 u 8:($3==11?($4==1?($5==1?($6==0.15?($7==1?$15:"1/0"):"1/0"):"1/0"):"1/0"):"1/0") title "{/Symbol m}=1.0" ls 1 lc rgb "blue"#, \
#	file1 u 8:($3==11?($4==1?($5==1?($6==0.15?($7==0.15?$15:"1/0"):"1/0"):"1/0"):"1/0"):"1/0") title "{/Symbol m}=0.15" ls 1 lc rgb "red"

#plot file1 u 6:($1==160211025057?($7==1?($8==1?$15:"1/0"):"1/0"):"1/0") title "{/Symbol m}=1.0" ls 1 lc rgb "blue" axes x1y1, \
#	file1 u 6:($1==160211120529?($7==1?($8==1?$15:"1/0"):"1/0"):"1/0") notitle ls 1 lc rgb "blue" axes x1y1, \
#	file1 u 6:($1==160211031204?($7==0.15?($8==1?$15:"1/0"):"1/0"):"1/0") title "{/Symbol m}=0.15" ls 1 lc rgb "red" axes x1y1, \
#	file1 u 6:($1==160211120434?($7==0.15?($8==1?$15:"1/0"):"1/0"):"1/0") notitle ls 1 lc rgb "red" axes x1y1

pause -1
