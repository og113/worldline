output="entropy_kappa_T_a_0.02_K_11_pot_8.dat"

reset
set xrange [0:10]
set yrange [0:2]
set isosamples 10, 10
set samples 10, 10
set cntrparam order 4
set cntrparam bspline
set dgrid3d
set table 'test.dat'
splot output u 1:2:3
unset table

set cont base
set cntrparam level incremental -4.0, 0.2, 4.0
unset surf
set table 'cont.dat'
splot output u 1:2:3
unset table

reset
set xrange [0:10]
set yrange [0:2]
set xlabel "{/Symbol k}"
set ylabel "T"
set zlabel "{/Symbol S}"
set title "{/Symbol S} vs {/Symbol k} and T: K=11, a=0.02"
unset key
set palette rgbformulae 33,13,10
l '<bash cont.sh cont.dat 0 7 0'
p 'test.dat' w ima, '<bash cont.sh cont.dat 1 7 0' w l lt -1 lw 1.5

pause -1
