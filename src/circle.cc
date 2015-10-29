/*
	circle
		program to generate circlular loop
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_trig.h>
#include <fstream>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - defining basic quantities
		2 - getting inputs from argv
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*-------------------------------------------------------------------------------------------------------------------------
	1 - defining basic quantities
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 4

ParametersRange pr;
pr.load("inputs2");
if (pr.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
Parameters p = pr.Min;
uint Length = pow(2,p.K);

/*-------------------------------------------------------------------------------------------------------------------------
	2 - getting inputs from argv
-------------------------------------------------------------------------------------------------------------------------*/

if (argc == 2) p.K = stringToNumber<uint>(argv[1]);
else if (argc % 2 && argc>1) {
for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("")==0);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

cout << "generating " << p.Nl << " circular unit loops each of " << Length << " points in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
string file, asciiFile;
uint Seed = time(NULL);
Loop<dim> loop(p.K,Seed);
Metropolis<dim> met(loop,p,Seed);
Point<dim> p0, point;
number R = 1.0, w = 2.0*pi/(number)Length;
if (abs(p.G)>MIN_NUMBER && abs(p.B)>MIN_NUMBER)
	R =  1.0/p.G/p.B;
	
for (uint j=0; j<p.Nl; j++) {
	file = "data/circle/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_R_"+nts<number>(R)\
														+"_run_"+nts<uint>(j)+".dat";
	for (uint k=0; k<Length; k++) {
		point = p0;
		point[2] += R*gsl_sf_cos(w*k);
		point[3] += R*gsl_sf_sin(w*k);
		loop[k] = point;
	}
	loop.save(file);
	loop.clear();
}

return 0;
}
