reset
set xrange [-5:0]
set yrange [2:5]
set isosample 150, 150
set table 'temp/test.dat'
splot sin(1.3*x)*cos(.9*y)+cos(.8*x)*sin(1.9*y)+cos(y*.2*x)
unset table

set cont base
set cntrparam level incremental -3, 0.5, 3
unset surf
set table 'temp/cont.dat'
splot sin(1.3*x)*cos(0.9*y)+cos(.8*x)*sin(1.9*y)+cos(y*.2*x)
unset table

reset
set xrange [-5:0]
set yrange [2:5]
unset key
set palette rgbformulae 33,13,10
l '<bash tests/contourLabel.sh temp/cont.dat 0 15 0'
p 'temp/test.dat' w ima, '<bash tests/contourLabel.sh temp/cont.dat 1 15 0' w l lt -1 lw 1.5

pause -1
