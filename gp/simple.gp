#gnuplot program to plot results from each group of loops, should give a gaussian plot

# setting print to stdout
set print "-"

# file to print from
if (exists("inFile")) {
	xLabel=system('echo "'.inFile.'" | sed -ne ''s/.*\/\([a-zA-Z]\+\)_\([a-zA-Z]\+\)\([0-9._ ]\+\).*/\2/p'' ');
	yLabel=system('echo "'.inFile.'" | sed -ne ''s/.*\/\([a-zA-Z]\+\)_\([a-zA-Z]\+\)\([0-9._ ]\+\).*/\1/p'' ');
	set title inFile noenhanced;
}
else {
	if (exists("inFiles")) {
		xLabel=system('echo "'.inFiles.'" | sed -ne ''s/.*\/\([a-zA-Z]\+\)_\([a-zA-Z]\+\)\([0-9._ ]\+\).*/\2/p'' ');
		yLabel=system('echo "'.inFiles.'" | sed -ne ''s/.*\/\([a-zA-Z]\+\)_\([a-zA-Z]\+\)\([0-9._ ]\+\).*/\1/p'' ');
		set title inFiles noenhanced;
	}
	else {
		print "must supply inFile or inFiles";
		exit;
	}
}

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

# setting title
if (exists("Title")) {
	set title Title noenhanced;
}

# setting logscale if required
if ((exists("logy")) || (yLabel eq "S")) \
set logscale y; \
set format y "10^{%L}"; \
else \
if (exists("logx")) \
set logscale x; \
set format x "10^{%L}"; \
else \
if (exists("log")) \
set logscale xy; \
set format xy "10^{%L}"; \

set key below
set autoscale
set xtic auto
set ytic auto
set xlabel xLabel
set ylabel yLabel

colours="blue red pink green black orange purple"

if (exists("inFile")) {
	if (exists("zero")) {
		plot inFile using 0:1 with lines ls 1 title inFile noenhanced;
	}
	else {
	plot inFile using 1:2 with lines ls 1 title inFile noenhanced;
	}
}
else {
	set for [i=1:words(colours)] lt 1 lc rgb word(colours, i);
	if (exists("zero")) {
		plot for [file in inFiles] file using 0:1 with lines title file noenhanced;
	}
	else {
		plot for [file in inFiles] file using 1:2 with lines title file noenhanced;
	}
}

if (!exists("outFile")) {
	pause -1
}
