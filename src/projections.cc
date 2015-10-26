/*
	projections
		prints all 2d projections of a loop for 4d
*/

#include <ctime>
#include <fstream>
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
		3 - loading and projecting loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*-------------------------------------------------------------------------------------------------------------------------
	1 - getting inputs from argv
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 4
string inputsFile = "inputs";
string loopFileIn = "data/s0/loops/dim_4/K_10/loop_run_0.dat";
string loopFileOut = "data/temp/projections.dat";

if (argc % 2 && argc>1) {
for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("fileIn")==0 || id.compare("fi")==0) loopFileIn = (string)argv[2*j+2];
		else if (id.compare("fileOut")==0 || id.compare("fo")==0) loopFileOut = (string)argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	2 - defining basic quantities
-------------------------------------------------------------------------------------------------------------------------*/

ParametersRange pr;
pr.load(inputsFile);
if (pr.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
Parameters p = pr.Min;
uint Length = pow(2,p.K);

cout << "forming 2d projections of loop from " << loopFileIn << endl;
cout << "saving to " << loopFileOut << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - loading and projecting loops
-------------------------------------------------------------------------------------------------------------------------*/

uint Seed = time(NULL);
Loop<dim> loop(p.K,Seed);
loop.load(loopFileIn);
ofstream os;
os.open(loopFileOut);
	
os << setprecision(16);
for (uint j=0; j<Length; j++) {
	for (uint k=0; k<dim; k++)
		os << setw(24) << (loop[j])[k];
	os << endl;
}

os.close();

return 0;
}
