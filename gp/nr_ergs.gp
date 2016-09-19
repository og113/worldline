# gnuplot program to plot ergs

file1="db/nr/plotsHand/ergs_a_K_11_kappa_1_pot_8_T_0.2.dat"
file2="db/nr/plotsHand/ergs_a_K_12_kappa_1_pot_8_T_0.2.dat"
file3="db/nr/plotsHand/ergs_a_K_11_kappa_0.1_pot_10_T_0.8.dat"
file5="db/nr/plotsHand/ergs_a_K_straight_11_kappa_1_pot_8_T_0.2.dat"
file6="db/nr/plotsHand/ergs_a_K_11_kappa_1_pot_9_T_0.2.dat"
file7="db/nr/plotsHand/ergs_a_K_12_kappa_1_pot_9_T_0.2.dat"

set key below
set autoscale
set xtic auto
set ytic auto
set title "Energy expressions versus a: {/Symbol k}=1, pot=8, T=0.2"
#set title "Energy expressions versus a: {/Symbol k}=0.1, pot=10, T=0.8"
set xlabel "a"
set ylabel "E"

# analytic approximate energy
pi=3.1415926535897932;
ergLowTemp(kappa,T)=2.0*(pi**4)*kappa*(T**5)/45.0
ergHighTemp(kappa)=2.0*(1.0-sqrt(kappa/4.0/pi))

#plot file5 using 1:5 with lines lt 1 lc rgb "blue" title "Local expression, K=11", \
#	file5 using 1:6 with lines lt 1 lc rgb "red" title "Thermal expression, K=11", \
#	file5 using 1:(ergHighTemp($2)) with lines lt 1 lc rgb "black" title "High temperature result", \

plot file1 using 1:5 with lines lt 1 lc rgb "blue" title "Local expression, K=11", \
	file1 using 1:6 with lines lt 1 lc rgb "red" title "Thermal expression, K=11", \
	file1 using 1:(ergLowTemp($2,$3)) with lines lt 1 lc rgb "black" title "Approximate low temperature result", \
	file2 using 1:5 with lines lt 1 lc rgb "orange" title "Local expression, K=12", \
	file2 using 1:6 with lines lt 1 lc rgb "green" title "Thermal expression, K=12", \
	file1 using 1:(-1.0/4.0/pi) with lines lt 1 lc rgb "black" title "1/4{/Symbol p}", \

pause -1

