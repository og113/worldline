#gnuplot program to plot worldlines

# setting print to stdout
set print "-"

# need to supply file/files to print
if (!exists("inFile") && !exists("inFiles")) \
print "must supply inFile or inFiles"; \
exit;

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

# stationary point of 3d potential
pi=3.1415926535897932
g=1.0
B=1.0
maximum(g,B)=sqrt(g/4.0/pi/B) # 1.771395468764214?
minimum=0.189673005548958 # for a=0.05

set key below
set autoscale
set xtic auto
set ytic auto

set title "Worldline, projection onto z-t plane"
set xlabel "z"
set ylabel "t"
set xrange [-1:1]
#set yrange [-1:1]

colours="blue red green orange cyan pink yellow black"

if (exists("max")) \
set arrow from (-maximum(g,B)/2.0),0 to (maximum(g,B)/2.0),0 nohead;

if (exists("min")) \
set arrow from (-minimum/2.0),1 to (minimum/2.0),1 nohead;

if (exists("inFile")) {
	plot inFile using 3:4 with lines ls 1 title inFile noenhanced;
}
else {
	set for [i=1:words(colours)] lt 1 lc rgb word(colours, i)
	plot for [file in inFiles] file using 3:4 with lines title file noenhanced;
}

if (!exists("outFile")) {
	pause -1
}
