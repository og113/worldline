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
		1 - getting inputs from argv
		2 - defining basic quantities
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string inputsFile = "inputs";
bool extend = false;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("extend")==0) extend = (stn<uint>(argv[2*j+2])!=0);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

cout << "using inputs file " << inputsFile << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	2 - defining basic quantities
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 4

ParametersRange pr;
pr.load(inputsFile);
if (pr.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
Parameters p = pr.Min;
uint Length = pow(2,p.K);

cout << "generating " << p.Nl << " circular unit loops each of " << Length << " points in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
string file, asciiFile;
uint Seed = time(NULL);
Loop<dim> loop(p.K,Seed);
Metropolis<dim> met(loop,p,Seed);
Point<dim> p0, point;
number R = 1.0, R0 = 1.0, w = 2.0*PI/(number)Length;
if (abs(p.G)>MIN_NUMBER && abs(p.B)>MIN_NUMBER) {
	R =  1.0/p.G/p.B;
	R0 = R;
}
if (abs(p.Epsi)>MIN_NUMBER && extend)
	R += p.Epsi;
	
for (uint j=0; j<p.Nl; j++) {
	file = "data/circle/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_R_"+nts<number>(R0)\
														+"_rank_"+nts<uint>(j)+".dat";
	if (j==0)
		cout << "printing to " << file << ", with runs from 0..." << p.Nl << endl;
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
