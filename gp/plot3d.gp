# gnuplot program to plot zmin

file="db/nr/plotsHand/zmin_B_T_K_11_pot_14.dat"

# if outFile exists, saving to file
if (exists("outFile")) {
	if (outFile ne 'gui') {
		set output outFile;
		suffix=system('echo "'.outFile.'" | sed -ne ''s/\(.*\)\.\(.*\)/\2/p'' ');
		if (suffix eq "png") {
			set term png size 1600,800;
		}
		else {
			if(suffix eq "svg") {
				set term svg size 1600,800;
			}
		}
	}
}

set key below
set autoscale
set xtic auto
set ytic auto
set y2tic auto
set title "z_{min}/a; a=0.01; pot=14"
set xlabel "{/Symbol k}"
set ylabel "T"

a=0.01


plot file using 1:2:($3/a) with points lt 1 palette

pause -1

