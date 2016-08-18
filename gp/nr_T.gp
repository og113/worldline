# gnuplot program to plot some early results of the disjoint temp calcs

file1="db/nr/plotsHand/S_T_kappa_0.01_K_11_pot_12.dat"
file2="db/nr/plotsHand/S_T_kappa_0.0316227766016838_K_11_pot_12.dat"
file3="db/nr/plotsHand/S_T_kappa_0.1_K_11_pot_12.dat"

set key below
set xtic auto
set ytic auto
#set logscale y
#set format y "10^{%L}"
set autoscale
set title "S vs T; a=0.003, K=11, pot=12"
set xlabel "T"
set ylabel "S"

# analytic quantities
pi=3.1415926535897932
Tmax(kappa)=sqrt(2.0)/(pi**0.75)/(kappa**0.25);
Sstraight(T,kappa)=2.0*(1.0-sqrt(kappa/4.0/pi))/T;

# values of kappa considered
kappa1=0.01
kappa2=0.0316227766016838
kappa3=0.1

# printing Tmax
#set arrow from (Tmax(kappa1)),2.0 to (Tmax(kappa1)),1.0 nohead

# titles
wavyTitle(kappa) = sprintf("{/Symbol k}=%g, wavy",kappa)
straightTitle(kappa) = sprintf("{/Symbol k}=%g, straight",kappa)

plot file1 using 1:2 with points title wavyTitle(kappa1) lt 1 lc rgb "red", \
	Sstraight(x,kappa1)with lines title straightTitle(kappa1) lt 1 lc rgb "red", \
	file2 using 1:2 with points title wavyTitle(kappa2) lt 1 lc rgb "blue", \
	Sstraight(x,kappa2) with lines title straightTitle(kappa2) lt 1 lc rgb "blue", \
	file3 using 1:2 with points title wavyTitle(kappa3) lt 1 lc rgb "green", \
	Sstraight(x,kappa3) with lines title straightTitle(kappa3) lt 1 lc rgb "green"

pause -1
