# gnuplot program to plot some early results of the disjoint temp calcs

file="db/nr/plotsHand/zmax_T_a_0.003_K_11_kappa_0.01_pot_12.dat"

set key below
set xtic auto
set ytic auto
#set logscale y
#set format y "10^{%L}"
set autoscale
set title "V3d vs T; a=0.003, K=11, pot=12"
set xlabel "T"
set ylabel "V3d"

# analytic quantities
pi=3.1415926535897932
V3d(kappa,a,r)=2.0-r-(kappa/4.0/pi)*(1.0/sqrt(r**2+a**2)-exp(-r**2/a**2)/a);
V3dStraight(kappa,a)=2.0*(1.0-sqrt(kappa/4.0/pi));

# values of parameters considered
kappa=0.01
a=0.003

plot file using 1:(V3d(kappa,a,$2)) with points title "V3d" lt 1 lc rgb "red", \
	V3dStraight(kappa,a) with lines title "straight line energy" lt 1 lc rgb "blue"

pause -1
