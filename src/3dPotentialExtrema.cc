/*-------------------------------------------------------------------------------------------------------------------------
	3dPotentialExtrema
		program to investigate the extrema of the 3d potential, as in email from Arttu on 22/3/16
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
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
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

using namespace std;

int main(int argc, char** argv) {
/*-------------------------------------------------------------------------------------------------------------------------
	1 - loading inputs and argv
-------------------------------------------------------------------------------------------------------------------------*/

// argv options
bool verbose = true;
bool bifurcate = true;
bool clever = false;
double frac = 0.5;
uint pot = 4;
string inputsFile = "inputs4";
string fo = "";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("bifurcate")==0) bifurcate = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("clever")==0) clever = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("frac")==0) frac = stn<double>(argv[2*j+2]);
		else if (id.compare("pot")==0) pot = stn<uint>(argv[2*j+2]);
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("out")==0 || id.compare("fo")==0) fo = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

// loading inputs
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;

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
uint Npl = (pr.Steps)[9];
if (Npl==0)
	Npl++;
uint minRuns = 2;
if (verbose && bifurcate) {
	cout << "a_min = " << a_min << ", a_max = " << a_max << endl;
	cout << "bisecting a maximum of " << Npl << " times" << endl;
	cout << "aiming for a tolerance of " << tol << endl;
}

// checking fo
if (fo.empty())
	fo = "data/nr/3dpot/extrema_bifurcate_"+nts((uint)bifurcate)+"_pot_"+nts(pot)+"_Kappa_"+nts(params.kappa)+".dat";

// initial guesses and ranges
double max_guess = sqrt(params.kappa/4.0/PI);
double max_l = max_guess - params.a;
double max_r = max_guess + params.a;
double min_guess = params.a;
double min_l = 0.0;
double min_r = max_l;

// asigning functions
gsl_function V_gsl;
V_gsl.params = &params;
gsl_function dV_gsl;
dV_gsl.params = &params;

if (pot==4) {
	V_gsl.function = &V4;
	dV_gsl.function = &dV4;
}
else if (pot==0) {
	V_gsl.function = &V0;
	dV_gsl.function = &dV0;
}
else {
	cerr << "3dPotentialExtrema Error: pot=" << pot << endl;
	return 1;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3 - finding minima and maxima
-------------------------------------------------------------------------------------------------------------------------*/

// objects to hold data
vector<double> a_vec(0), max_vec(0), min_vec(0), Vmax_vec(0), Vmin_vec(0);
double Vmin_l = 0.0, Vmin_r = 0.0, test = 1.0;

// results file
ofstream os;
os.open(fo.c_str());
if (os.good()) {
	os << left << setprecision(16);
}
else {
	cerr << "3dPotentialExtrema error: cannot write to " << fo << endl;
	os.close();
	return 1;
}

if (verbose) {
	printf("\n");
	printf("%6s%6s%6s%20s%20s%20s%20s%20s%20s\n","run","pot","kappa","a","max approx","max","min","Vmax","Vmin");
}

uint run = 0;
while((test>tol || run<minRuns) && run<Npl) {
	
	// stepping parameters
	if (bifurcate) {
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
	}
	else {
		if (run==0) {
			params.a = a_min;
		}
		else {
			p = pr.position(run);
			params.setFromParameters(p);
		}
		if (clever) {
			max_l = ((max_guess - 3.0*pow(params.a,2))>min_guess? (max_guess - 3.0*pow(params.a,2)):\
							 0.5*((1.0-frac)*max_guess+frac*min_guess) );
			max_r = max_guess + 3.0*pow(params.a,2);
		}
		else {
			max_l = max_guess - params.a;
			max_r = max_guess + params.a;
		}
		min_l = 0.0;
		min_r = max_l;
	}
	
	// setting a_vec
	a_vec.push_back(params.a);	
	
	// finding extrema
	max_vec.push_back(brentRootFinder(&dV_gsl,max_guess,max_l,max_r));
	min_vec.push_back(brentRootFinder(&dV_gsl,min_guess,min_l,min_r));

	// finding value of V at min and max
	Vmax_vec.push_back(V_gsl.function(max_vec[run],&params));
	Vmin_vec.push_back(V_gsl.function(min_vec[run],&params));
	
	// test
	test = (bifurcate? abs(Vmin_vec[run]): 1.0);
	
	if (bifurcate) {
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
	}
	else {
		max_guess = max_vec[run];
		min_guess = min_vec[run];
	}
	
	os << setw(25) << pot << setw(25) << params.kappa << setw(25) << params.a;
	os << setw(25) << max_vec[run] << setw(25) << min_vec[run];
	os << setw(25) << Vmax_vec[run] << setw(25) << Vmin_vec[run] << endl;
	
	number rGuess = sqrt(params.kappa/4.0/PI) - 3.0*sqrt(PI/4.0/params.kappa)*pow(params.a,2)\
						 - 15.0*pow(PI/params.kappa,3.0/2.0)*pow(params.a,4)/4.0;
	if (verbose)
		printf("%6i%6i%6.2g%20.12g%20.12g%20.12g%20.12g%20.12g%20.12g\n",run,pot,params.kappa,params.a,\
			rGuess,max_vec[run],min_vec[run],Vmax_vec[run],Vmin_vec[run]);
		
	run++;
}

if (os.good())
	os.close();

if (verbose)
	printf("\n");

cout << "printed results: " << fo << endl;

return 0;
}
