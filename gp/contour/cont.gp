reset
set xrange [0.:1]
set yrange [0:2]
#set xrange [0.01:0.15]
#set yrange [0.8:1.9]
set isosamples 10, 10
set samples 10, 10
set cntrparam order 4
set cntrparam bspline
set dgrid3d
set table 'test.dat'
splot 'output.dat' u 1:2:3
unset table

set cont base
set cntrparam level incremental 0, 0.1, 3.2
unset surf
set table 'cont.dat'
splot 'output.dat' u 1:2:3
unset table

reset
set xrange [0.:1]
set yrange [0:2]
#set xrange [0.01:0.15]
#set yrange [0.8:1.9]
unset key
set palette rgbformulae 33,13,10
l '<bash cont.sh cont.dat 0 10 0'
p 'test.dat' w ima, '<bash cont.sh cont.dat 1 10 0' w l lt -1 lw 1.5

pause -1
