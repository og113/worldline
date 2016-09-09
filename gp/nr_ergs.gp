# gnuplot program to plot ergs

file="db/nr/plotsHand/ergs_a_K_11_kappa_1_pot_8_T_0.2.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Energy expressions versus a: K=11, {/Symbol k}=1, pot=8, T=0.2"
set xlabel "a"
set ylabel "E"

# analytic approximate energy
pi=3.1415926535897932;
ergLowTemp(kappa,T)=2.0*(pi**4)*kappa*(T**5)/45.0
ergHighTemp(kappa)=2.0*(1.0-sqrt(kappa/4.0/pi))

plot file using 1:5 with lines lt 1 lc rgb "blue" title "Local expression", \
	file using 1:6 with lines lt 1 lc rgb "red" title "Thermal expression", \
	file using 1:(ergLowTemp($2,$3)) with lines lt 1 lc rgb "green" title "Approximate low temperature result", \



pause -1

