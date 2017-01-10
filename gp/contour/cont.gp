output8="/home/oli/Physics/worldline/db/nr/plotsHand/sExtrapolated_kappa_T_K_12_pot_13.dat"
output10="/home/oli/Physics/worldline/db/nr/plotsHand/sExtrapolated_kappa_T_K_12_pot_15.dat"

pi=3.1415926535897932
Tmax(kappa)=sqrt(2)/pi**0.75/kappa**0.25
Sstraight(kappa,T)=(1.0-sqrt(kappa/4.0/pi))*2.0/T
Tstraight(kappa,s)=(1.0-sqrt(kappa/4.0/pi))*2.0/s
min2(a,b)=(a<b?a:b);
min3(a,b,c)=min2(a,min2(b,c));

reset
set xrange [0:10]
set yrange [0:2]
set isosamples 10, 10
set samples 10, 10
set cntrparam order 4
set cntrparam bspline
set dgrid3d
set table '/home/oli/Physics/worldline/gp/contour/test8.dat'
splot output8 u 1:2:3
unset table

set cont base
set cntrparam level incremental 0.0, 0.1, 4.0
unset surf
set table '/home/oli/Physics/worldline/gp/contour/cont8.dat'
splot output8 u 1:2:3
unset table

reset
set xrange [0:10]
set yrange [0:2]
set isosamples 10, 10
set samples 10, 10
set cntrparam order 4
set cntrparam bspline
set dgrid3d
set table '/home/oli/Physics/worldline/gp/contour/test10.dat'
splot output10 u 1:2:3
unset table

set cont base
set cntrparam level incremental 0.0, 0.1, 4.0
unset surf
set table '/home/oli/Physics/worldline/gp/contour/cont10.dat'
splot output10 u 1:2:3
unset table

reset
set xrange [0:1]
set yrange [0:2]
set xlabel "{/Symbol k}"
set ylabel "T"
set zlabel "S"
set title "S vs {/Symbol k} and T: K=12"
set key below
set palette rgbformulae 33,13,10
l '<bash /home/oli/Physics/worldline/gp/contour/cont.sh /home/oli/Physics/worldline/gp/contour/cont8.dat 0 2 1'
l '<bash /home/oli/Physics/worldline/gp/contour/cont.sh /home/oli/Physics/worldline/gp/contour/cont10.dat 0 2 1'

#p 'test.dat' w ima, \
#	'<bash /home/oli/Physics/worldline/gp/contour/cont.sh /home/oli/Physics/worldline/gp/contour/cont.dat 1 2 1' w l lt -1 lw 1.5, \
#	Tmax(x) w lines dt 2 lt 1 lc rgb "black" title "{/Symbol S}=0"
	
p '<bash /home/oli/Physics/worldline/gp/contour/cont.sh /home/oli/Physics/worldline/gp/contour/cont8.dat 1 2 1' w l lt -1 lw 1.5 lc rgb "black", \
	'<bash /home/oli/Physics/worldline/gp/contour/cont.sh /home/oli/Physics/worldline/gp/contour/cont10.dat 1 2 1' w l lt -1 lw 1.5 lc rgb "red", \
	Tmax(x) w lines dt 2 lt 1 lc rgb "green" title "{/Symbol S}=0", \
	Tstraight(x,1.2) w lines lt 1 lc rgb "blue" title "S=1.2", \
	Tstraight(x,1.4) w lines lt 1 lc rgb "blue" title "S=1.4", \
	Tstraight(x,1.6) w lines lt 1 lc rgb "orange" title "S=1.6", \
	Tstraight(x,1.8) w lines lt 1 lc rgb "cyan" title "S=1.8", \
	Tstraight(x,2.0) w lines lt 1 lc rgb "magenta" title "S=2.0", \
	Tstraight(x,2.2) w lines lt 1 lc rgb "grey" title "S=2.2"

pause -1
