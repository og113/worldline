/*-------------------------------------------------------------------------------------------------------------------------
	3dPotentialExtrema
		program to investigate the extrema of the 3d potential, as in email from Arttu on 22/3/16
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "simple.h"
#include "gsl_extras.h"
#include "3dpotential.h"
#include "print.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - loading inputs and argv
	2 - setting up
	3 - parameter loop, finding minima and maxima
	4 - printing results
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/


using namespace std;

int main(int argc, char** argv) {
/*-------------------------------------------------------------------------------------------------------------------------
	1 - loading inputs and argv
-------------------------------------------------------------------------------------------------------------------------*/

// argv options
bool verbose = true;
string inputsFile = "inputs4";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

// loading inputs
ParametersRange pr;
pr.load(inputsFile);

if (verbose) {
	cout << "3dPotentialExtrema:" << endl << endl;
	cout << "using inputs file " << inputsFile << ":" << endl;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2 - setting up
-------------------------------------------------------------------------------------------------------------------------*/

// parameters
params_for_V params;
params.setFromParameters(pr.Min);
double a_min = (pr.Min).Epsi;
double a_max = (pr.Max).Epsi;
double tol = (pr.Min).Lambda;
uint Npl = pr.totalSteps();
uint minRuns = 2;
if (verbose) {
	cout << "a_min = " << a_min << ", a_max = " << a_max << endl;
	cout << "bisecting a maximum of " << Npl << " times" << endl;
	cout << "aiming for a tolerance of " << tol << endl;
}

// initial guesses and ranges
double max_guess = sqrt(params.kappa/4.0/PI);
double max_l = max_guess - params.a;
double max_r = max_guess + params.a;
double min_guess = params.a;
double min_l = 0.0;
double min_r = max_l;

// asigning functions
gsl_function dV_gsl;
dV_gsl.function = &dV;
dV_gsl.params = &params;

/*-------------------------------------------------------------------------------------------------------------------------
	2 - finding minima and maxima
-------------------------------------------------------------------------------------------------------------------------*/

// objects to hold data
vector<double> a_vec(0), max_vec(0), min_vec(0), Vmax_vec(0), Vmin_vec(0);
double Vmin_l = 0.0, Vmin_r = 0.0, test = 1.0;

if (verbose) {
	printf("\n");
	printf("%6s%6s%24s%24s%24s%24s%24s\n","run","kappa","a","max","min","Vmax","Vmin");
}

uint run = 0;

while((test>tol || run<minRuns) && run<Npl) {
	
	// stepping parameters
	if (run==0) {
		params.a = a_min;
	}
	else if (run==1) {
		params.a = a_max;
	}
	else if (run>1) {
		params.a = 0.5*(a_min + a_max);
	}
	
	max_l = max_guess - params.a;
	max_r = max_guess + params.a;
	min_guess = params.a;
	min_l = 0.0;
	min_r = max_l;
	
	// setting a_vec
	a_vec.push_back(params.a);	
	
	// finding extrema
	max_vec.push_back(brentRootFinder(&dV_gsl,max_guess,max_l,max_r));
	min_vec.push_back(brentRootFinder(&dV_gsl,min_guess,min_l,min_r));

	// finding value of V at min and max
	Vmax_vec.push_back(V(max_vec[run],&params));
	Vmin_vec.push_back(V(min_vec[run],&params));
	
	// test
	test = abs(Vmin_vec[run]);
	
	if (run==0) {
		Vmin_l = Vmin_vec[0];
	}
	else if (run==1) {
		Vmin_r = Vmin_vec[1];
		// checking to see if change of sign over range of 'a'
		if (Vmin_l*Vmin_r>0.0) {
			cerr << "error: Vmin_l*Vmin_r>0.0)" << endl;
			cerr << "change range of 'a' to search" << endl;
			break;
			return 1;
		}
	}
	else if (run>1) {
		// doing bifurcation
		if (Vmin_vec[run]*Vmin_l<0.0) {
			a_max = params.a;
		}
		else {
			a_min = params.a;
		}
	}
	
	if (verbose)
		printf("%6i%6.2g%24.16g%24.16g%24.16g%24.16g%24.16g\n",run,params.kappa,params.a,\
			max_vec[run],min_vec[run],Vmax_vec[run],Vmin_vec[run]);
		
	run++;
}

if (verbose)
	printf("\n");

/*-------------------------------------------------------------------------------------------------------------------------
	3 - printing results	
-------------------------------------------------------------------------------------------------------------------------*/

string fo = "data/nr/3dpot/extrema_Kappa_"+nts(params.kappa)+".dat";
saveVectorAscii(fo,a_vec);
saveVectorAsciiAppend(fo,max_vec);
saveVectorAsciiAppend(fo,min_vec);
saveVectorAsciiAppend(fo,Vmax_vec);
saveVectorAsciiAppend(fo,Vmin_vec);


return 0;
}
